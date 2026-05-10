"""Python translation of analysis/classes/MERFISHDecoder.m.

This class coordinates MERFISH data normalization, preprocessing, decoding,
segmentation, molecule/barcode access, feature quantification, reports, and
persistence. The Python implementation follows the MATLAB public method layout
and uses numpy/pandas/scikit-image/tifffile equivalents for the portable logic.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import json
import math
import pickle
import re
import shutil

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import exposure, filters, measure, morphology, segmentation, transform, restoration
import tifffile

from .FoundFeature import FoundFeature


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    suffix = path.suffix.lower()
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path)
    if suffix in {".tsv"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path)


def _write_table(df: pd.DataFrame, path: Path) -> None:
    _ensure_dir(path.parent)
    if path.suffix.lower() == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)


def _gaussian_kernel(size: int, sigma: float) -> np.ndarray:
    ax = np.arange(-size // 2 + 1, size // 2 + 1)
    xx, yy = np.meshgrid(ax, ax)
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    return kernel / np.sum(kernel)


def _parse_barcode_string(value: Any) -> List[int]:
    if isinstance(value, (list, tuple, np.ndarray)):
        return [int(v) for v in value]
    text = str(value).strip()
    if not text:
        return []
    nums = re.findall(r"[01]", text)
    return [int(x) for x in nums]


@dataclass
class MERFISHDecoder:
    rawDataPath: str | Path
    normalizedDataPath: str | Path
    verbose: bool = True
    overwrite: bool = False
    sliceIDs: list = field(default_factory=list)
    name: str = ""
    hal_version: str = "hal1"
    dataOrganizationPath: str | Path | None = None
    codebookPath: str | Path | None = None
    parameters: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    version: str = "0.6"
    mDecoderPath: str = "mDecoder"
    fiducialDataPath: str = "fiducial_data"
    warpedDataPath: str = "warped_data"
    processedDataPath: str = "processed_data"
    barcodePath: str = "barcodes"
    reportPath: str = "reports"
    segmentationPath: str = "segmentation"
    summationPath: str = "summation"
    mosaicPath: str = "mosaics"

    dataOrganization: pd.DataFrame = field(default_factory=pd.DataFrame)
    codebook: pd.DataFrame = field(default_factory=pd.DataFrame)
    codebookHeader: Dict[str, Any] = field(default_factory=dict)
    mappedRawData: pd.DataFrame = field(default_factory=pd.DataFrame)
    scaleFactors: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    fovIDs: List[int] = field(default_factory=list)
    fovPos: np.ndarray = field(default_factory=lambda: np.empty((0, 2), dtype=float))
    zPos: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    imageSize: Tuple[int, int] = (0, 0)
    bitNames: List[str] = field(default_factory=list)
    numBits: int = 0
    numBarcodes: int = 0
    numFov: int = 0
    pixelSize: float = 109.0
    parallel: Any = None
    numPar: int = 1

    def __post_init__(self) -> None:
        self.rawDataPath = Path(self.rawDataPath)
        self.normalizedDataPath = Path(self.normalizedDataPath)
        for sub in [
            self.mDecoderPath,
            self.fiducialDataPath,
            self.warpedDataPath,
            self.processedDataPath,
            self.barcodePath,
            self.reportPath,
            self.segmentationPath,
            self.summationPath,
            self.mosaicPath,
        ]:
            _ensure_dir(self.normalizedDataPath / sub)
        if not self.parameters:
            self.parameters, _ = self.InitializeParameters()
        if self.dataOrganizationPath:
            self.dataOrganization = _read_table(Path(self.dataOrganizationPath))
        if self.codebookPath:
            self._load_codebook(Path(self.codebookPath))
        self.LoadImageMetaData()
        if self.scaleFactors.size == 0 and self.numBits:
            self.scaleFactors = np.ones(self.numBits, dtype=float)

    def _load_codebook(self, path: Path) -> None:
        self.codebookPath = path
        rows = []
        header = {}
        if path.exists():
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                if not line.strip():
                    continue
                if line.startswith("#"):
                    parts = line[1:].split(":", 1)
                    if len(parts) == 2:
                        header[parts[0].strip()] = parts[1].strip()
                    continue
                rows.append(line)
        if rows:
            from io import StringIO
            text = "\n".join(rows)
            try:
                df = pd.read_csv(StringIO(text))
            except Exception:
                df = pd.read_csv(StringIO(text), sep="\t")
        else:
            df = pd.DataFrame()
        self.codebook = df
        self.codebookHeader = header
        if not df.empty:
            bit_cols = [c for c in df.columns if re.fullmatch(r"bit\d+|Bit\d+|b\d+|B\d+", str(c))]
            if not bit_cols and "barcode" in df.columns:
                max_len = max(len(_parse_barcode_string(v)) for v in df["barcode"])
                for i in range(max_len):
                    df[f"bit{i+1}"] = [(_parse_barcode_string(v) + [0] * max_len)[i] for v in df["barcode"]]
                bit_cols = [f"bit{i+1}" for i in range(max_len)]
            self.bitNames = bit_cols
            self.numBits = len(bit_cols)
            self.numBarcodes = len(df)

    def _subpath(self, name: str) -> Path:
        return self.normalizedDataPath / name

    def FoundFeaturesToCSV(self, outputPath: str | Path | None = None, **kwargs) -> pd.DataFrame:
        features = self.GetFoundFeatures()
        frames = [f.Feature2Table() for f in features]
        table = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        path = Path(outputPath) if outputPath else self._subpath(self.segmentationPath) / "found_features.csv"
        _write_table(table, path)
        return table

    def BarcodesToCSV(self, outputPath: str | Path | None = None, parsed: bool = False, **kwargs) -> pd.DataFrame:
        rows = []
        for fov in self.fovIDs or self._discover_fovs():
            df = self.GetParsedBarcodeList(fov) if parsed else self.GetBarcodeList(fov)
            if not df.empty:
                rows.append(df.assign(fovID=fov))
        out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        path = Path(outputPath) if outputPath else self._subpath(self.barcodePath) / ("parsed_barcodes.csv" if parsed else "barcodes.csv")
        _write_table(out, path)
        return out

    def CalculateFeatureCounts(self, outputPath: str | Path | None = None, **kwargs) -> pd.DataFrame:
        features = self.GetFoundFeatures()
        barcode_df = self.BarcodesToCSV(parsed=True)
        if barcode_df.empty or not features:
            counts = pd.DataFrame()
        else:
            if "feature_uID" in barcode_df.columns:
                group_cols = ["feature_uID", "barcode_id" if "barcode_id" in barcode_df.columns else "barcodeID"]
                counts = barcode_df.groupby(group_cols).size().reset_index(name="count")
            else:
                rows = []
                for f in features:
                    for _, r in barcode_df.iterrows():
                        z = r.get("z", r.get("zIndex", 0))
                        if f.IsInFeature([r.get("x", r.get("abs_x", np.nan)), r.get("y", r.get("abs_y", np.nan)), z]):
                            rows.append({"feature_uID": f.uID, "barcode_id": r.get("barcode_id", r.get("barcodeID", np.nan)), "count": 1})
                counts = pd.DataFrame(rows).groupby(["feature_uID", "barcode_id"]).size().reset_index(name="count") if rows else pd.DataFrame()
        path = Path(outputPath) if outputPath else self._subpath(self.segmentationPath) / "feature_counts.csv"
        _write_table(counts, path)
        return counts

    def GenerateFoundFeaturesReport(self) -> pd.DataFrame:
        df = self.FoundFeaturesToCSV()
        report = df.describe(include="all") if not df.empty else pd.DataFrame()
        _write_table(report.reset_index(), self._subpath(self.reportPath) / "found_features_report.csv")
        return report

    def GenerateSummationReport(self) -> pd.DataFrame:
        data = self.GetSummedSignal()
        normalizedSignal, sumSignal, sumPixels, channelNames, featureUIDs = data
        df = pd.DataFrame(sumSignal, columns=channelNames if channelNames else None)
        if len(featureUIDs) == len(df):
            df.insert(0, "feature_uID", featureUIDs)
        _write_table(df, self._subpath(self.reportPath) / "summation_report.csv")
        return df

    def GenerateFeatureCountsReport(self) -> pd.DataFrame:
        counts = self.CalculateFeatureCounts()
        if counts.empty:
            report = counts
        else:
            report = counts.groupby("feature_uID")["count"].sum().reset_index(name="total_count")
        _write_table(report, self._subpath(self.reportPath) / "feature_counts_report.csv")
        return report

    def LoadImageMetaData(self) -> None:
        if not self.dataOrganization.empty:
            cols = {c.lower(): c for c in self.dataOrganization.columns}
            if "bitname" in cols:
                self.bitNames = list(pd.unique(self.dataOrganization[cols["bitname"]].astype(str)))
            elif "name" in cols:
                self.bitNames = list(pd.unique(self.dataOrganization[cols["name"]].astype(str)))
            self.numBits = max(self.numBits, len(self.bitNames))
        self.mappedRawData = self.MapRawData()
        if not self.mappedRawData.empty:
            if "fovID" in self.mappedRawData.columns:
                self.fovIDs = sorted(int(x) for x in pd.unique(self.mappedRawData["fovID"].dropna()))
            self.numFov = len(self.fovIDs)
            first = Path(self.mappedRawData.iloc[0]["path"])
            try:
                arr = tifffile.imread(first)
                if arr.ndim >= 2:
                    self.imageSize = tuple(int(v) for v in arr.shape[-2:])
            except Exception:
                return None

    def MapRawData(self) -> pd.DataFrame:
        rows = []
        if not self.rawDataPath.exists():
            self.mappedRawData = pd.DataFrame()
            return self.mappedRawData
        pattern = re.compile(r"(?:fov|FOV|f)(\d+)")
        for p in sorted(self.rawDataPath.rglob("*")):
            if p.suffix.lower() in {".tif", ".tiff", ".dax", ".npy"}:
                m = pattern.search(p.name)
                fov = int(m.group(1)) if m else len(rows) + 1
                rows.append({"path": str(p), "file": p.name, "fovID": fov, "extension": p.suffix.lower()})
        self.mappedRawData = pd.DataFrame(rows)
        return self.mappedRawData

    def SegmentFOV(self, fovIDs: Iterable[int] | int) -> List[FoundFeature]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        all_features = []
        params = self.parameters.get("segmentation", {})
        min_size = int(params.get("minCellSize", 100))
        for fov in fovs:
            img, _ = self.GetImage(int(fov), dataChannel=None, zStack=None)
            if img.size == 0:
                continue
            frame = img.astype(float)
            if frame.ndim == 3:
                frame2 = np.max(frame, axis=0) if frame.shape[0] < min(frame.shape[-2:]) else np.max(frame, axis=2)
            else:
                frame2 = frame
            smooth = filters.gaussian(frame2, sigma=max(float(params.get("watershedFrameFilterSize", 5)), 0) / 3.0)
            thresh = filters.threshold_otsu(smooth) if np.any(smooth) else 0
            mask = smooth > thresh
            mask = morphology.remove_small_objects(mask, min_size=min_size)
            labels = measure.label(mask)
            features = []
            for lab in range(1, int(labels.max()) + 1):
                ff = FoundFeature(labels, fov, self.GetFOVCoordinates(fov), self.pixelSize, params.get("stageOrientation", [1, 1]), params.get("boundingBox", [-100, -100, 200, 200]), [0], lab)
                ff.AssignFeatureID(lab)
                features.append(ff)
            all_features.extend(features)
            with open(self._subpath(self.segmentationPath) / f"found_features_fov_{int(fov)}.pkl", "wb") as fh:
                pickle.dump(features, fh)
        return all_features

    def GenerateLowResolutionMosaic(self) -> np.ndarray:
        fovs = self.fovIDs or self._discover_fovs()
        if not fovs:
            return np.empty((0, 0))
        down = int(self.parameters.get("display", {}).get("downSample", 10))
        images, coords = [], []
        for fov in fovs:
            img, xy = self.GetImage(fov, None, None)
            if img.size == 0:
                continue
            frame = np.squeeze(img)
            while frame.ndim > 2:
                frame = frame.max(axis=0)
            small = transform.resize(frame, (max(1, frame.shape[0] // down), max(1, frame.shape[1] // down)), anti_aliasing=True)
            images.append(small)
            coords.append(xy)
        if not images:
            return np.empty((0, 0))
        h, w = images[0].shape
        n = int(math.ceil(math.sqrt(len(images))))
        mosaic = np.zeros((n * h, n * w), dtype=float)
        for i, im in enumerate(images):
            r, c = divmod(i, n)
            mosaic[r*h:(r+1)*h, c*w:(c+1)*w] = im
        tifffile.imwrite(self._subpath(self.mosaicPath) / "low_res_mosaic.tif", mosaic.astype(np.float32))
        return mosaic

    def SumRawSignalFOV(self, fovIDs: Iterable[int] | int):
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        outputs = []
        for fov in fovs:
            img, _ = self.GetImage(int(fov), None, None)
            if img.size == 0:
                continue
            summed = np.sum(img, axis=0) if img.ndim > 2 else img
            out = self._subpath(self.summationPath) / f"raw_sum_fov_{int(fov)}.npy"
            np.save(out, summed)
            outputs.append(out)
        return outputs

    def CombineRawSum(self) -> np.ndarray:
        arrays = [np.load(p) for p in sorted(self._subpath(self.summationPath).glob("raw_sum_fov_*.npy"))]
        combined = np.sum(arrays, axis=0) if arrays else np.array([])
        if combined.size:
            np.save(self._subpath(self.summationPath) / "raw_sum_combined.npy", combined)
        return combined

    def GenerateCompositeImage(self):
        return self.GenerateLowResolutionMosaic()

    def GetFOVCoordinates(self, fovID: int) -> np.ndarray:
        if self.fovPos.size and self.fovIDs and int(fovID) in self.fovIDs:
            return self.fovPos[self.fovIDs.index(int(fovID))]
        if not self.dataOrganization.empty:
            cols = {c.lower(): c for c in self.dataOrganization.columns}
            if "fovid" in cols and {"x", "y"}.issubset(cols):
                row = self.dataOrganization[self.dataOrganization[cols["fovid"]] == int(fovID)]
                if not row.empty:
                    return row[[cols["x"], cols["y"]]].iloc[0].to_numpy(float)
        return np.array([0.0, 0.0], dtype=float)

    def GetMosaic(self, sliceID: int | None = None, framesToLoad: Sequence[int] | None = None):
        path = self._subpath(self.mosaicPath) / "low_res_mosaic.tif"
        if not path.exists():
            image = self.GenerateLowResolutionMosaic()
        else:
            image = tifffile.imread(path)
        return image, np.array([self.GetFOVCoordinates(f) for f in (self.fovIDs or [])])

    def GetSummedSignal(self):
        path = self._subpath(self.summationPath) / "summed_signal.csv"
        if path.exists():
            df = pd.read_csv(path)
            channel_names = [c for c in df.columns if c != "feature_uID"]
            values = df[channel_names].to_numpy(float) if channel_names else np.empty((len(df), 0))
            norm = values / np.maximum(values.sum(axis=1, keepdims=True), 1)
            return norm, values, values.sum(axis=1), channel_names, list(df.get("feature_uID", []))
        return np.empty((0, 0)), np.empty((0, 0)), np.array([]), [], []

    def ParseFOV(self, fovIDs: Iterable[int] | int) -> List[pd.DataFrame]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        features = self.GetFoundFeatures()
        outputs = []
        for fov in fovs:
            barcodes = self.GetBarcodeList(int(fov))
            if barcodes.empty:
                outputs.append(barcodes)
                continue
            rows = []
            for _, r in barcodes.iterrows():
                x = r.get("abs_x", r.get("x", np.nan))
                y = r.get("abs_y", r.get("y", np.nan))
                z = r.get("z", 0)
                assigned = ""
                for feat in features:
                    if feat.InFov([fov]) and feat.IsInFeature([x, y, z]):
                        assigned = feat.uID
                        break
                rr = r.to_dict()
                rr["feature_uID"] = assigned
                rows.append(rr)
            out = pd.DataFrame(rows)
            _write_table(out, self._subpath(self.barcodePath) / f"parsed_barcodes_fov_{int(fov)}.csv")
            outputs.append(out)
        return outputs

    def FindMoleculesFOV(self, fovIDs: Iterable[int] | int) -> List[pd.DataFrame]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        outputs = []
        p = self.parameters.get("molecules", {})
        thresh = float(p.get("molIntensityThreshold", 1000))
        sigma = float(p.get("molLowPassfilterSize", 5)) / 3.0
        for fov in fovs:
            img, _ = self.GetImage(int(fov), None, None)
            frame = np.squeeze(img)
            while frame.ndim > 2:
                frame = frame.max(axis=0)
            filt = filters.gaussian(frame.astype(float), sigma=sigma)
            labels = measure.label(filt > thresh)
            props = measure.regionprops_table(labels, intensity_image=frame, properties=("label", "centroid", "area", "max_intensity", "mean_intensity"))
            df = pd.DataFrame(props).rename(columns={"centroid-0": "y", "centroid-1": "x"})
            _write_table(df, self._subpath(self.processedDataPath) / f"molecules_fov_{int(fov)}.csv")
            outputs.append(df)
        return outputs

    def CombineFeatures(self) -> List[FoundFeature]:
        features = self.GetFoundFeatures()
        out = self._subpath(self.segmentationPath) / "found_features_combined.pkl"
        with open(out, "wb") as fh:
            pickle.dump(features, fh)
        return features

    def GetFoundFeatures(self) -> List[FoundFeature]:
        combined = self._subpath(self.segmentationPath) / "found_features_combined.pkl"
        if combined.exists():
            with open(combined, "rb") as fh:
                return pickle.load(fh)
        features = []
        for p in sorted(self._subpath(self.segmentationPath).glob("found_features_fov_*.pkl")):
            with open(p, "rb") as fh:
                features.extend(pickle.load(fh))
        return features

    def _discover_fovs(self) -> List[int]:
        if self.fovIDs:
            return self.fovIDs
        if self.mappedRawData.empty:
            self.MapRawData()
        return sorted(int(x) for x in pd.unique(self.mappedRawData.get("fovID", pd.Series(dtype=int)).dropna()))

    def _raw_file_for_fov(self, fovID: int) -> Optional[Path]:
        if self.mappedRawData.empty:
            self.MapRawData()
        rows = self.mappedRawData[self.mappedRawData["fovID"] == int(fovID)] if not self.mappedRawData.empty else pd.DataFrame()
        if rows.empty:
            return None
        return Path(rows.iloc[0]["path"])

    def GetImage(self, fovID: int, dataChannel: Any = None, zStack: Any = None):
        p = self._raw_file_for_fov(int(fovID))
        if p is None or not p.exists():
            return np.array([]), self.GetFOVCoordinates(fovID)
        if p.suffix.lower() == ".npy":
            img = np.load(p)
        else:
            img = tifffile.imread(p)
        if dataChannel is not None and img.ndim >= 3:
            img = img[int(dataChannel)]
        if zStack is not None and img.ndim >= 3:
            img = img[..., int(zStack)]
        return img, self.GetFOVCoordinates(fovID)

    def GetProcessedImage(self, fovID: int, dataChannel: Any = None, zStack: Any = None):
        path = self._subpath(self.processedDataPath) / f"processed_fov_{int(fovID)}.tif"
        if path.exists():
            img = tifffile.imread(path)
            return img, self.GetFOVCoordinates(fovID)
        return self.GetImage(fovID, dataChannel, zStack)

    def GetDecodedImage(self, fovID: int):
        d = self._subpath(self.barcodePath) / f"decoded_fov_{int(fovID)}.tif"
        m = self._subpath(self.barcodePath) / f"magnitude_fov_{int(fovID)}.tif"
        decoded = tifffile.imread(d) if d.exists() else np.array([])
        mag = tifffile.imread(m) if m.exists() else np.array([])
        return decoded, mag, self.GetFOVCoordinates(fovID)

    def GetBarcodeList(self, fovID: int) -> pd.DataFrame:
        path = self._subpath(self.barcodePath) / f"barcodes_fov_{int(fovID)}.csv"
        return pd.read_csv(path) if path.exists() else pd.DataFrame()

    def GetMoleculeList(self, fovID: int) -> pd.DataFrame:
        path = self._subpath(self.processedDataPath) / f"molecules_fov_{int(fovID)}.csv"
        return pd.read_csv(path) if path.exists() else pd.DataFrame()

    def GetParsedBarcodeList(self, fovID: int) -> pd.DataFrame:
        path = self._subpath(self.barcodePath) / f"parsed_barcodes_fov_{int(fovID)}.csv"
        return pd.read_csv(path) if path.exists() else pd.DataFrame()

    def WarpFOV(self, fovIDs: Iterable[int] | int) -> List[Path]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        out_paths = []
        for fov in fovs:
            img, _ = self.GetImage(int(fov), None, None)
            if img.size == 0:
                continue
            out = self._subpath(self.warpedDataPath) / f"warped_fov_{int(fov)}.tif"
            tifffile.imwrite(out, img)
            out_paths.append(out)
        return out_paths

    def PreprocessFOV(self, fovIDs: Iterable[int] | int) -> List[Path]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        params = self.parameters.get("preprocess", {})
        method = params.get("preprocessingMethod", "highPassDecon")
        out_paths = []
        for fov in fovs:
            img, _ = self.GetImage(int(fov), None, None)
            if img.size == 0:
                continue
            arr = img.astype(float)
            sigma = float(params.get("highPassKernelSize", 3))
            low = ndi.gaussian_filter(arr, sigma=sigma)
            processed = arr - low
            processed -= processed.min()
            if "decon" in str(method).lower():
                kernel = _gaussian_kernel(10, 2)
                if processed.ndim == 2:
                    processed = restoration.richardson_lucy(processed, kernel, num_iter=int(params.get("numIterDecon", 20)))
            out = self._subpath(self.processedDataPath) / f"processed_fov_{int(fov)}.tif"
            tifffile.imwrite(out, processed.astype(np.float32))
            out_paths.append(out)
        return out_paths

    def CheckStatus(self, fovIDs: Iterable[int] | int, analysisType: str) -> bool:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        paths = {
            "warp": self.warpedDataPath,
            "preprocess": self.processedDataPath,
            "decode": self.barcodePath,
            "segment": self.segmentationPath,
            "sum": self.summationPath,
        }
        folder = self._subpath(paths.get(str(analysisType).lower(), str(analysisType)))
        return all(any(folder.glob(f"*{int(fov)}*")) for fov in fovs)

    def GenerateWarpReport(self) -> pd.DataFrame:
        files = list(self._subpath(self.warpedDataPath).glob("*.tif"))
        df = pd.DataFrame({"file": [str(p) for p in files], "size_bytes": [p.stat().st_size for p in files]})
        _write_table(df, self._subpath(self.reportPath) / "warp_report.csv")
        return df

    def ResetScaleFactors(self) -> None:
        self.scaleFactors = np.ones(self.numBits or 1, dtype=float)

    def OptimizeScaleFactors(self, numFOV: int | None = None, **kwargs) -> np.ndarray:
        fovs = (self.fovIDs or self._discover_fovs())[: int(numFOV or self.parameters.get("optimization", {}).get("optNumFov", 50))]
        traces = []
        for fov in fovs:
            img, _ = self.GetProcessedImage(fov, None, None)
            if img.size == 0:
                img, _ = self.GetImage(fov, None, None)
            arr = np.asarray(img, dtype=float)
            if arr.ndim >= 3:
                flat = arr.reshape((-1, arr.shape[-1])) if arr.shape[-1] == self.numBits else arr.reshape((arr.shape[0], -1)).T
                traces.append(flat)
        if traces:
            data = np.vstack(traces)
            q = float(self.parameters.get("optimization", {}).get("quantileTarget", 0.9))
            sf = np.nanquantile(data, q, axis=0)
            sf[sf <= 0] = 1.0
            self.scaleFactors = sf
        else:
            self.ResetScaleFactors()
        return self.scaleFactors

    def SetScaleFactors(self, scaleFactors: Sequence[float]) -> None:
        sf = np.asarray(scaleFactors, dtype=float)
        if sf.ndim != 1:
            raise ValueError("scaleFactors must be a 1D vector")
        self.scaleFactors = sf
        self.numBits = len(sf)

    def DecodeFOV(self, fovIDs: Iterable[int] | int) -> List[pd.DataFrame]:
        fovs = [fovIDs] if isinstance(fovIDs, (int, np.integer)) else list(fovIDs)
        vectors, single = self.GenerateDecodingMatrices()
        outputs = []
        for fov in fovs:
            stack = self.ReadAndFilterTiffStack(self._raw_file_for_fov(int(fov)))
            if stack.size == 0:
                outputs.append(pd.DataFrame())
                continue
            decoded, mag, traces, D = self.DecodePixels(stack, self.scaleFactors, vectors, self.parameters.get("decoding", {}).get("distanceThreshold", 0.5176))
            self.SaveDecodedImageAndMagnitudeImage(decoded, mag, int(fov))
            df = self.GenerateBarcodes({}, mag, single, traces, D, None, int(fov), decodedImage=decoded)
            _write_table(df, self._subpath(self.barcodePath) / f"barcodes_fov_{int(fov)}.csv")
            outputs.append(df)
        return outputs

    def SaveDecodedImageAndMagnitudeImage(self, decodedImage: np.ndarray, localMagnitude: np.ndarray, fovID: int) -> None:
        _ensure_dir(self._subpath(self.barcodePath))
        tifffile.imwrite(self._subpath(self.barcodePath) / f"decoded_fov_{int(fovID)}.tif", np.asarray(decodedImage).astype(np.uint16))
        mag = np.asarray(localMagnitude)
        if decodedImage.size and mag.ndim == 1:
            mag = mag.reshape(np.asarray(decodedImage).shape)
        tifffile.imwrite(self._subpath(self.barcodePath) / f"magnitude_fov_{int(fovID)}.tif", mag.astype(np.float32))

    def Pixel2Abs(self, pixelPos: Sequence[float] | np.ndarray, fovID: int) -> np.ndarray:
        pts = np.asarray(pixelPos, dtype=float)
        one = pts.ndim == 1
        pts = pts.reshape(1, -1) if one else pts
        center = self.GetFOVCoordinates(fovID)
        pixel_um = float(self.parameters.get("decoding", {}).get("pixelSize", self.pixelSize)) / 1000.0
        shape = np.array(self.imageSize if self.imageSize != (0, 0) else [0, 0], dtype=float)
        x = (pts[:, 0] - shape[1] / 2.0) * pixel_um + center[0]
        y = (pts[:, 1] - shape[0] / 2.0) * pixel_um + center[1]
        out = np.column_stack([x, y])
        return out[0] if one else out

    def GenerateBarcodes(self, properties, localMagnitude, singleBitErrorBarcodes, pixelTraces, D, b, fovID, decodedImage=None) -> pd.DataFrame:
        if decodedImage is None:
            return pd.DataFrame()
        labels = np.asarray(decodedImage)
        rows = []
        for barcode_id in sorted(int(x) for x in np.unique(labels) if int(x) > 0):
            comp = measure.label(labels == barcode_id)
            props = measure.regionprops(comp, intensity_image=np.asarray(localMagnitude).reshape(labels.shape) if np.asarray(localMagnitude).size == labels.size else None)
            for p in props:
                y, x = p.centroid[:2]
                abs_xy = self.Pixel2Abs([x, y], fovID)
                rows.append({
                    "barcode_id": barcode_id,
                    "area": p.area,
                    "x": x,
                    "y": y,
                    "abs_x": abs_xy[0],
                    "abs_y": abs_xy[1],
                    "fovID": fovID,
                    "mean_intensity": getattr(p, "mean_intensity", np.nan),
                    "min_distance": float(np.nanmin(D)) if np.asarray(D).size else np.nan,
                })
        return pd.DataFrame(rows)

    def ReadAndFilterTiffStack(self, tiffName2Read) -> np.ndarray:
        if tiffName2Read is None:
            return np.array([])
        path = Path(tiffName2Read)
        if not path.exists():
            return np.array([])
        arr = np.load(path) if path.suffix.lower() == ".npy" else tifffile.imread(path)
        arr = np.asarray(arr, dtype=float)
        low = float(self.parameters.get("decoding", {}).get("lowPassKernelSize", 1))
        if low > 0:
            arr = ndi.uniform_filter(arr, size=tuple([1] * max(arr.ndim - 2, 0) + [int(low), int(low)])) if arr.ndim >= 2 else arr
        crop = int(self.parameters.get("decoding", {}).get("crop", 0))
        if crop > 0 and arr.ndim >= 2 and min(arr.shape[-2:]) > 2 * crop:
            arr = arr[..., crop:-crop, crop:-crop]
        return arr

    def InitializeScaleFactors(self) -> np.ndarray:
        self.scaleFactors = np.ones(self.numBits or 1, dtype=float)
        return self.scaleFactors

    def GenerateDecodingMatrices(self):
        if self.codebook.empty:
            if self.numBits == 0:
                self.numBits = len(self.scaleFactors) if self.scaleFactors.size else 1
            eye = np.eye(self.numBits, dtype=float)
            return eye, eye.copy()
        bit_cols = [c for c in self.codebook.columns if re.fullmatch(r"bit\d+|Bit\d+|b\d+|B\d+", str(c))]
        if not bit_cols and "barcode" in self.codebook.columns:
            bits = [_parse_barcode_string(v) for v in self.codebook["barcode"]]
            max_len = max(len(v) for v in bits)
            mat = np.array([v + [0] * (max_len - len(v)) for v in bits], dtype=float)
        else:
            mat = self.codebook[bit_cols].to_numpy(float)
        norms = np.linalg.norm(mat, axis=1)
        norms[norms == 0] = 1
        vectors = mat / norms[:, None]
        error_words = []
        for row in mat.astype(int):
            for i in range(row.size):
                rr = row.copy()
                rr[i] = 1 - rr[i]
                error_words.append(rr)
        single = np.unique(np.asarray(error_words, dtype=float), axis=0) if error_words else np.empty((0, mat.shape[1]))
        n2 = np.linalg.norm(single, axis=1) if single.size else np.array([])
        if single.size:
            n2[n2 == 0] = 1
            single = single / n2[:, None]
        self.numBits = mat.shape[1]
        self.numBarcodes = mat.shape[0]
        if self.scaleFactors.size != self.numBits:
            self.scaleFactors = np.ones(self.numBits, dtype=float)
        return vectors, single

    def GenerateOptimizationReport(self, localScaleFactors=None, onBitIntensity=None, allCounts=None) -> pd.DataFrame:
        df = pd.DataFrame({"scale_factor": np.asarray(localScaleFactors if localScaleFactors is not None else self.scaleFactors, dtype=float)})
        _write_table(df, self._subpath(self.reportPath) / "optimization_report.csv")
        return df

    def SetParallel(self, p: Any) -> None:
        self.parallel = p
        self.numPar = getattr(p, "_processes", 1) if p is not None else 1

    def Save(self, path: str | Path | None = None, **kwargs) -> Path:
        out = Path(path) if path else self._subpath(self.mDecoderPath) / "MERFISHDecoder.pkl"
        _ensure_dir(out.parent)
        with open(out, "wb") as fh:
            pickle.dump(self, fh)
        return out

    def SetParameter(self, **kwargs) -> None:
        for key, value in kwargs.items():
            found = False
            for group in self.parameters.values():
                if key in group:
                    group[key] = value
                    found = True
            if not found:
                self.parameters.setdefault("custom", {})[key] = value

    def UpdateNormalizedDataPath(self, newPath: str | Path) -> None:
        self.normalizedDataPath = Path(newPath)
        for sub in [self.mDecoderPath, self.fiducialDataPath, self.warpedDataPath, self.processedDataPath, self.barcodePath, self.reportPath, self.segmentationPath, self.summationPath, self.mosaicPath]:
            _ensure_dir(self.normalizedDataPath / sub)

    def UpdateField(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)

    def LoadField(self, fieldName: str):
        return getattr(self, fieldName)

    def Downsample(self, bitNamesToRemove: Sequence[str] = (), fovIDsToRemove: Sequence[int] = ()) -> None:
        remove_bits = set(bitNamesToRemove)
        if remove_bits and self.bitNames:
            keep = [b for b in self.bitNames if b not in remove_bits]
            self.bitNames = keep
            self.numBits = len(keep)
        remove_fovs = set(int(x) for x in fovIDsToRemove)
        if remove_fovs:
            self.fovIDs = [f for f in self.fovIDs if f not in remove_fovs]
            self.numFov = len(self.fovIDs)

    @staticmethod
    def Load(dirPath: str | Path, **kwargs) -> "MERFISHDecoder":
        path = Path(dirPath)
        if path.is_dir():
            path = path / "mDecoder" / "MERFISHDecoder.pkl"
        with open(path, "rb") as fh:
            return pickle.load(fh)

    @staticmethod
    def DefaultWarpParameters():
        return {
            "warpingDataPath": None,
            "fiducialFitMethod": "daoSTORM",
            "controlPointOffsetRange": np.arange(-60, 60.5, 0.5),
            "numNN": 10,
            "pairDistanceTolerance": 3,
            "pixelSize": 109,
            "sigmaInit": 1.6,
            "daoThreshold": 500,
            "daoBaseline": 100,
            "exportWarpedBeads": True,
            "cameraOrientation": [0, 0, 0],
            "geoTransformEdges": [np.arange(0, 2048 + 25, 25), np.arange(0, 2048 + 25, 25)],
            "colorTransforms": [],
        }

    @staticmethod
    def DefaultPreprocessingParameters():
        return {
            "preprocessingMethod": "highPassDecon",
            "highPassKernelSize": 3,
            "deconKernel": _gaussian_kernel(10, 2),
            "numIterDecon": 20,
            "erosionElement": morphology.disk(1),
        }

    @staticmethod
    def DefaultDecodingParameters():
        return {
            "lowPassKernelSize": 1,
            "crop": 40,
            "decodingMethod": "distanceHS1",
            "distanceThreshold": 0.5176,
            "minBrightness": 1.0,
            "minArea": 1,
            "connectivity": np.ones((3, 3, 3), dtype=bool),
            "stageOrientation": [1, 1],
            "pixelSize": 109,
        }

    @staticmethod
    def DefaultSegmentationParameters():
        return {
            "segmentationMethod": "seededWatershed",
            "watershedSeedChannel": "DAPI",
            "watershedChannel": "polyT",
            "seedFrameFilterSize": 5,
            "seedThreshold": "adaptive",
            "minCellSize": 100,
            "watershedFrameFilterSize": 5,
            "watershedFrameThreshold": "adaptive",
            "ignoreZ": False,
            "boundingBox": [-100, -100, 200, 200],
            "maxEdgeDistance": 4,
            "maxFeatureCentroidDistance": 5,
            "dilationSize": 0.1,
            "saveSegmentationReports": False,
        }

    @staticmethod
    def DefaultSummationParameters():
        return {"sumExactOnly": False, "sumNormalizeByPixelCount": True}

    @staticmethod
    def DefaultMoleculeParameters():
        return {
            "molLowPassfilterSize": 5,
            "molIntensityThreshold": 1000,
            "molNumPixelSum": 1,
            "molDataChannels": [["RS0763", 4], ["RS1199", 4], ["RS1040", 4]],
        }

    @staticmethod
    def DefaultOptimizationParameters():
        return {
            "weightingOptimizationMethod": "equalOnBits",
            "quantileTarget": 0.9,
            "areaThresh": 4,
            "optNumFov": 50,
            "numIterOpt": 10,
            "normalizeToOne": False,
        }

    @staticmethod
    def DefaultDisplayParameters():
        return {"visibleOption": "on", "overwrite": True, "formats": ["fig", "png"], "useExportFig": False, "downSample": 10, "mosaicZInd": 4}

    @staticmethod
    def DefaultQuantificationParameters():
        return {"minimumBarcodeArea": 4, "minimumBarcodeBrightness": 10**0.75, "minimumDistanceToFeature": float("inf"), "zSliceRange": []}

    @staticmethod
    def InitializeParameters(**kwargs):
        defaults = {
            "warp": MERFISHDecoder.DefaultWarpParameters(),
            "preprocess": MERFISHDecoder.DefaultPreprocessingParameters(),
            "decoding": MERFISHDecoder.DefaultDecodingParameters(),
            "optimization": MERFISHDecoder.DefaultOptimizationParameters(),
            "display": MERFISHDecoder.DefaultDisplayParameters(),
            "segmentation": MERFISHDecoder.DefaultSegmentationParameters(),
            "quantification": MERFISHDecoder.DefaultQuantificationParameters(),
            "summation": MERFISHDecoder.DefaultSummationParameters(),
            "molecules": MERFISHDecoder.DefaultMoleculeParameters(),
        }
        for key, value in kwargs.items():
            placed = False
            for group in defaults.values():
                if key in group:
                    group[key] = value
                    placed = True
            if not placed:
                defaults.setdefault("custom", {})[key] = value
        return defaults, defaults.copy()

    @staticmethod
    def CheckTiffStack(tiffFileName: str | Path, expectedNumFrames: int) -> bool:
        try:
            with tifffile.TiffFile(tiffFileName) as tif:
                return len(tif.pages) < int(expectedNumFrames)
        except Exception:
            return True

    @staticmethod
    def DecodePixels(imageStack: np.ndarray, scaleFactors: Sequence[float], decodingVectors: np.ndarray, distanceThreshold: float):
        arr = np.asarray(imageStack, dtype=float)
        if arr.ndim == 3:
            imageHeight, imageWidth, stackLength = arr.shape
            numZPos = None
            pixelTraces = arr.reshape((imageHeight * imageWidth, stackLength))
            out_shape = (imageHeight, imageWidth)
        elif arr.ndim == 4:
            imageHeight, imageWidth, numZPos, stackLength = arr.shape
            pixelTraces = arr.reshape((imageHeight * imageWidth * numZPos, stackLength))
            out_shape = (imageHeight, imageWidth, numZPos)
        else:
            raise ValueError("imageStack must have shape (H,W,bits) or (H,W,Z,bits)")
        sf = np.asarray(scaleFactors, dtype=float)
        if sf.size != pixelTraces.shape[1]:
            sf = np.ones(pixelTraces.shape[1], dtype=float)
        sf[sf == 0] = 1.0
        pixelTraces = pixelTraces / sf[None, :]
        localMagnitude = np.sqrt(np.sum(pixelTraces * pixelTraces, axis=1))
        good = localMagnitude > 0
        pixelTraces[good] = pixelTraces[good] / localMagnitude[good, None]
        vectors = np.asarray(decodingVectors, dtype=float)
        if vectors.ndim != 2:
            raise ValueError("decodingVectors must be 2D")
        tree = cKDTree(vectors)
        D, barcodeID = tree.query(pixelTraces, k=1)
        exact = D <= float(distanceThreshold)
        decoded = np.zeros(pixelTraces.shape[0], dtype=np.uint16)
        decoded[exact] = barcodeID[exact] + 1
        decodedImage = decoded.reshape(out_shape)
        return decodedImage, localMagnitude.reshape(out_shape), pixelTraces, D
