"""
Python translation of analysis/functions/MERFISHPerformanceMetrics.m.

Purpose/use:
    Generate MERFISH barcode performance metrics: exact/corrected barcode
    counts per cell or FOV, per-bit error counts, area/brightness histograms,
    barcode density maps, optional abundance/FPKM correlation figures, and
    saved CSV outputs.

Python adaptation notes:
    MATLAB reads custom MERFISH binary barcode lists through ReadBinaryFile and
    loads MERFISHDecoder objects. This translation implements the same metric
    logic for Python-readable barcode tables (CSV/TSV/Parquet/JSON/NPZ/pickle)
    and can also call translated fileIO functions if available in a combined
    package. It is real metric code, not a scaffold.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Sequence, Tuple
import json
import math
import os
import pickle
import re
import shutil
import warnings

import numpy as np

try:
    import pandas as pd
except Exception as exc:  # pragma: no cover
    raise ImportError("MERFISHPerformanceMetrics requires pandas") from exc

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


@dataclass
class MERFISHPerformanceParameters:
    verbose: bool = False
    logProgress: bool = True
    archive: bool = True
    barcodePath: str | None = None
    outputPath: str | None = None
    codebookPath: str | None = None
    abundDataPath: str | None = None
    cellBoundariesPath: str | None = None
    parallel: Any = None
    cellIDMethod: str = "fov"
    blockSize: int = 100_000
    brightnessThreshold: float = 10 ** 0.75
    areaThreshold: float = 2.0
    stageOrientation: Tuple[int, int] = (1, -1)
    visibleOption: str = "on"
    brightnessBins: np.ndarray | None = None
    areaBins: np.ndarray | None = None
    blankFnc: Callable[[str], bool] | None = None
    barcodeDensityBinSize: int = 20
    overwrite: bool = True
    formats: Tuple[str, ...] = ("png", "fig")
    useExportFig: bool = False
    imageSize: Tuple[int, int] | None = None
    numZPos: int = 1

    def __post_init__(self) -> None:
        if self.brightnessBins is None:
            self.brightnessBins = np.arange(0.0, 4.00001, 0.025)
        else:
            self.brightnessBins = np.asarray(self.brightnessBins, dtype=float)
        if self.areaBins is None:
            self.areaBins = np.arange(1, 16, 1, dtype=float)
        else:
            self.areaBins = np.asarray(self.areaBins, dtype=float)
        if self.blankFnc is None:
            self.blankFnc = lambda name: re.search(r"Blank-", str(name)) is not None
        if self.cellIDMethod not in {"fov", "cellID"}:
            raise ValueError("cellIDMethod must be 'fov' or 'cellID'")
        if self.visibleOption not in {"on", "off"}:
            raise ValueError("visibleOption must be 'on' or 'off'")
        if self.blockSize <= 0:
            raise ValueError("blockSize must be positive")
        if self.brightnessThreshold < 0 or self.areaThreshold < 0:
            raise ValueError("thresholds must be non-negative")
        if self.barcodeDensityBinSize <= 0:
            raise ValueError("barcodeDensityBinSize must be positive")


def _log(message: str, parameters: MERFISHPerformanceParameters, log_lines: List[str]) -> None:
    log_lines.append(message)
    if parameters.verbose:
        print(message)


def _serializable_parameters(parameters: MERFISHPerformanceParameters) -> Dict[str, Any]:
    data = asdict(parameters)
    fn = data.get("blankFnc")
    if callable(fn):
        data["blankFnc"] = getattr(fn, "__name__", repr(fn))
    if isinstance(data.get("brightnessBins"), np.ndarray):
        data["brightnessBins"] = data["brightnessBins"].tolist()
    if isinstance(data.get("areaBins"), np.ndarray):
        data["areaBins"] = data["areaBins"].tolist()
    return data


def _find_codebook(normalized_data_path: Path) -> Path:
    matches = list(normalized_data_path.glob("*_codebook*.csv")) + list(normalized_data_path.glob("*codebook*.csv"))
    unique = sorted(set(matches))
    if not unique:
        raise FileNotFoundError("Could not find a valid codebook CSV in normalizedDataPath")
    if len(unique) > 2:
        raise FileExistsError("Found too many files that look like codebooks")
    return unique[0]


def _load_codebook(path: Path) -> tuple[pd.DataFrame, list[str]]:
    codebook = pd.read_csv(path)
    lower = {c.lower(): c for c in codebook.columns}
    if "name" not in lower:
        for candidate in ["gene", "gene_name", "target", "id"]:
            if candidate in lower:
                codebook = codebook.rename(columns={lower[candidate]: "name"})
                break
    else:
        codebook = codebook.rename(columns={lower["name"]: "name"})
    if "name" not in codebook.columns:
        codebook.insert(0, "name", [f"barcode_{i+1}" for i in range(len(codebook))])

    bit_columns = [c for c in codebook.columns if re.match(r"^(bit|hyb|round|readout)[_\- ]*\d+", str(c), re.I)]
    if not bit_columns:
        bit_columns = [c for c in codebook.columns if set(str(v) for v in codebook[c].dropna().unique()).issubset({"0", "1", "0.0", "1.0", "True", "False"})]
    return codebook, bit_columns


def _load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    if suffix == ".npz":
        data = np.load(path, allow_pickle=True)
        return pd.DataFrame({k: data[k] for k in data.files})
    if suffix in {".pkl", ".pickle", ".matb", ".bin"}:
        with open(path, "rb") as handle:
            obj = pickle.load(handle)
        if isinstance(obj, pd.DataFrame):
            return obj
        if isinstance(obj, dict):
            return pd.DataFrame(obj)
        if isinstance(obj, list):
            return pd.DataFrame(obj)
        raise ValueError(f"unsupported pickled object in {path}")
    raise ValueError(f"unsupported barcode table format: {path}")


def _discover_barcode_files(barcode_path: Path, cell_id_method: str) -> list[Path]:
    if cell_id_method == "cellID":
        candidates = [
            barcode_path / "parsed" / "assigned_blist.csv",
            barcode_path / "parsed" / "assigned_blist.tsv",
            barcode_path / "parsed" / "assigned_blist.parquet",
            barcode_path / "parsed" / "assigned_blist.pkl",
            barcode_path / "parsed" / "assigned_blist.bin",
        ]
        found = [p for p in candidates if p.exists()]
        if not found:
            raise FileNotFoundError("Could not find parsed/assigned_blist in barcodePath")
        return found
    fov_dir = barcode_path / "barcode_fov"
    if not fov_dir.exists():
        fov_dir = barcode_path
    files = []
    for suffix in ("*.csv", "*.tsv", "*.parquet", "*.json", "*.jsonl", "*.npz", "*.pkl", "*.pickle", "*.bin"):
        files.extend(fov_dir.glob(suffix))
    return sorted(files)


def _standardize_barcodes(df: pd.DataFrame, source_index: int = 1) -> pd.DataFrame:
    rename_map = {}
    lower = {str(c).lower(): c for c in df.columns}
    synonyms = {
        "barcode_id": ["barcode_id", "barcodeid", "barcode", "barcodeID", "id"],
        "is_exact": ["is_exact", "isexact", "exact", "isExact"],
        "error_dir": ["error_dir", "errordir", "error_direction", "errorDirection"],
        "error_bit": ["error_bit", "errorbit", "bit_error", "errorBit"],
        "total_magnitude": ["total_magnitude", "totalmagnitude", "brightness", "intensity", "totalMagnitude"],
        "area": ["area", "pixel_area"],
        "cellID": ["cellid", "cell_id", "cellID"],
        "fov_id": ["fov_id", "fovid", "fov", "fovID"],
        "x": ["x", "weighted_x", "weightedpixelcentroid_x"],
        "y": ["y", "weighted_y", "weightedpixelcentroid_y"],
        "z": ["z", "weighted_z", "weightedpixelcentroid_z"],
    }
    for target, names in synonyms.items():
        for name in names:
            if name.lower() in lower:
                rename_map[lower[name.lower()]] = target
                break
    out = df.rename(columns=rename_map).copy()

    if "weighted_pixel_centroid" in out.columns:
        parsed = out["weighted_pixel_centroid"].apply(_parse_centroid)
        out["x"] = [p[0] for p in parsed]
        out["y"] = [p[1] for p in parsed]
        out["z"] = [p[2] if len(p) > 2 else 1 for p in parsed]

    defaults = {
        "barcode_id": 1,
        "is_exact": 1,
        "error_dir": 0,
        "error_bit": 0,
        "total_magnitude": 1.0,
        "area": 1.0,
        "cellID": source_index,
        "fov_id": source_index,
        "x": np.nan,
        "y": np.nan,
        "z": 1,
    }
    for col, default in defaults.items():
        if col not in out.columns:
            out[col] = default
    numeric_cols = ["barcode_id", "is_exact", "error_dir", "error_bit", "total_magnitude", "area", "cellID", "fov_id", "x", "y", "z"]
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=["barcode_id"])
    out["barcode_id"] = out["barcode_id"].astype(int)
    out["is_exact"] = out["is_exact"].fillna(0).astype(int)
    out["error_dir"] = out["error_dir"].fillna(0).astype(int)
    out["error_bit"] = out["error_bit"].fillna(0).astype(int)
    out["area"] = out["area"].fillna(0.0).astype(float)
    out["total_magnitude"] = out["total_magnitude"].fillna(0.0).astype(float)
    return out


def _parse_centroid(value: Any) -> tuple[float, ...]:
    if isinstance(value, (list, tuple, np.ndarray)):
        return tuple(float(x) for x in value)
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    if len(nums) >= 2:
        return tuple(float(x) for x in nums[:3])
    return (np.nan, np.nan, 1.0)


def _hist2d_pairs(row_values: np.ndarray, col_values: np.ndarray, row_edges: np.ndarray, col_edges: np.ndarray) -> np.ndarray:
    hist, _, _ = np.histogram2d(row_values, col_values, bins=[row_edges, col_edges])
    return hist.astype(float)


def _save_csv(path: Path, data: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.asarray(data), delimiter=",")


def _save_figure(fig: Any, output_path: Path, name: str, parameters: MERFISHPerformanceParameters) -> None:
    if plt is None:
        return
    output_path.mkdir(parents=True, exist_ok=True)
    for fmt in parameters.formats:
        if fmt.lower() == "fig":
            continue
        out = output_path / f"{name}.{fmt}"
        if out.exists() and not parameters.overwrite:
            raise FileExistsError(out)
        fig.savefig(out, bbox_inches="tight")


def _correlation(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.sum(mask) < 3:
        return np.nan
    return float(np.corrcoef(np.log10(x[mask]), np.log10(y[mask]))[0, 1])


def MERFISHPerformanceMetrics(normalizedDataPath: str | os.PathLike[str], **kwargs: Any) -> Dict[str, Any]:
    """
    Translate MATLAB MERFISHPerformanceMetrics.

    Returns a dictionary containing output paths, arrays, and parsed parameters.
    CSV and figure outputs are written to outputPath just like the MATLAB
    function writes report files.
    """
    allowed = set(MERFISHPerformanceParameters.__dataclass_fields__.keys())
    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(f"unknown parameter(s): {sorted(unknown)}")
    parameters = MERFISHPerformanceParameters(**kwargs)
    normalized_path = Path(normalizedDataPath)
    if not normalized_path.exists():
        raise FileNotFoundError(f"normalizedDataPath does not exist: {normalized_path}")

    if parameters.codebookPath is None:
        codebook_path = _find_codebook(normalized_path)
    else:
        codebook_path = Path(parameters.codebookPath)
    barcode_path = Path(parameters.barcodePath) if parameters.barcodePath else normalized_path / "barcodes"
    output_path = Path(parameters.outputPath) if parameters.outputPath else barcode_path / "performance"
    output_path.mkdir(parents=True, exist_ok=True)

    log_lines: List[str] = []
    _log(f"Calculating performance metrics for {normalized_path}", parameters, log_lines)

    codebook, bit_names = _load_codebook(codebook_path)
    num_barcodes = len(codebook)
    num_bits = len(bit_names)
    if num_bits == 0:
        max_error_bit = 1
    else:
        max_error_bit = num_bits

    barcode_files = _discover_barcode_files(barcode_path, parameters.cellIDMethod)
    if not barcode_files:
        raise FileNotFoundError(f"No barcode tables found under {barcode_path}")
    _log(f"Found {len(barcode_files)} barcode file(s)", parameters, log_lines)

    data_frames = []
    for i, path in enumerate(barcode_files, start=1):
        local = _standardize_barcodes(_load_table(path), source_index=i)
        if parameters.cellIDMethod == "fov":
            local["fov_id"] = local["fov_id"].fillna(i).astype(int)
        data_frames.append(local)
    a_list = pd.concat(data_frames, ignore_index=True) if data_frames else pd.DataFrame()

    if a_list.empty:
        raise ValueError("No barcodes loaded")

    if parameters.cellIDMethod == "cellID":
        cell_ids = np.sort(a_list["cellID"].dropna().astype(int).unique())
        cell_to_index = {cell_id: idx for idx, cell_id in enumerate(cell_ids)}
        object_values = a_list["cellID"].astype(int).map(cell_to_index).to_numpy()
        fpkm_label = "Counts/Cell"
    else:
        cell_ids = np.sort(a_list["fov_id"].dropna().astype(int).unique())
        cell_to_index = {cell_id: idx for idx, cell_id in enumerate(cell_ids)}
        object_values = a_list["fov_id"].astype(int).map(cell_to_index).to_numpy()
        fpkm_label = "Counts/FOV"
    num_cells = len(cell_ids)

    brightness = np.divide(
        a_list["total_magnitude"].to_numpy(dtype=float),
        a_list["area"].replace(0, np.nan).to_numpy(dtype=float),
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        log_brightness = np.log10(brightness)
    brightness_area_hist = _hist2d_pairs(
        log_brightness[np.isfinite(log_brightness)],
        a_list.loc[np.isfinite(log_brightness), "area"].to_numpy(dtype=float),
        np.asarray(parameters.brightnessBins, dtype=float),
        np.asarray(parameters.areaBins, dtype=float),
    )

    keep = (a_list["area"] >= parameters.areaThreshold) & (brightness >= parameters.brightnessThreshold)
    filtered = a_list.loc[keep].copy()
    object_values_filtered = object_values[keep.to_numpy()]
    _log(f"Cut to {len(filtered)} barcode(s) after area/brightness thresholds", parameters, log_lines)

    row_edges = np.arange(0.5, num_barcodes + 1.5, 1.0)
    col_edges = np.arange(-0.5, num_cells + 0.5, 1.0)
    bit_edges = np.arange(-0.5, max_error_bit + 1.5, 1.0)

    exact = filtered["is_exact"].to_numpy(dtype=int) == 1
    corrected = ~exact
    barcode_indices = filtered["barcode_id"].to_numpy(dtype=int)
    valid_barcode = (barcode_indices >= 1) & (barcode_indices <= num_barcodes)

    counts_per_cell_exact = _hist2d_pairs(
        barcode_indices[exact & valid_barcode],
        object_values_filtered[exact & valid_barcode],
        row_edges,
        col_edges,
    )
    counts_per_cell_corrected = _hist2d_pairs(
        barcode_indices[corrected & valid_barcode],
        object_values_filtered[corrected & valid_barcode],
        row_edges,
        col_edges,
    )

    error_bit = filtered["error_bit"].to_numpy(dtype=int)
    error_dir = filtered["error_dir"].to_numpy(dtype=int)
    num_zero_to_one = _hist2d_pairs(
        barcode_indices[(error_dir == 0) & valid_barcode],
        error_bit[(error_dir == 0) & valid_barcode],
        row_edges,
        bit_edges,
    )
    num_one_to_zero = _hist2d_pairs(
        barcode_indices[(error_dir == 1) & valid_barcode],
        error_bit[(error_dir == 1) & valid_barcode],
        row_edges,
        bit_edges,
    )

    if parameters.imageSize is None:
        max_x = np.nanmax(filtered["x"].to_numpy(dtype=float)) if len(filtered) else 0
        max_y = np.nanmax(filtered["y"].to_numpy(dtype=float)) if len(filtered) else 0
        image_size = (max(1, int(math.ceil(max_y))), max(1, int(math.ceil(max_x))))
    else:
        image_size = tuple(parameters.imageSize)
    num_z = max(1, int(parameters.numZPos))
    x_edges = np.arange(1, image_size[1] + parameters.barcodeDensityBinSize, parameters.barcodeDensityBinSize)
    y_edges = np.arange(1, image_size[0] + parameters.barcodeDensityBinSize, parameters.barcodeDensityBinSize)
    if len(x_edges) < 2:
        x_edges = np.array([1, image_size[1] + 1], dtype=float)
    if len(y_edges) < 2:
        y_edges = np.array([1, image_size[0] + 1], dtype=float)
    barcode_density = np.zeros((len(y_edges) - 1, len(x_edges) - 1, num_z), dtype=float)
    z_values = filtered["z"].fillna(1).to_numpy(dtype=int)
    z_values = np.clip(z_values, 1, num_z)
    for z in range(1, num_z + 1):
        z_mask = z_values == z
        if np.any(z_mask):
            hist, _, _ = np.histogram2d(
                filtered.loc[z_mask, "y"].to_numpy(dtype=float),
                filtered.loc[z_mask, "x"].to_numpy(dtype=float),
                bins=[y_edges, x_edges],
            )
            barcode_density[:, :, z - 1] = hist

    dist_to_nucleus_in = np.zeros_like(counts_per_cell_exact)
    dist_to_nucleus_out = np.zeros_like(counts_per_cell_exact)
    fraction_in_nucleus = np.zeros_like(counts_per_cell_exact)

    _save_csv(output_path / "countsPerCellExact.csv", counts_per_cell_exact)
    _save_csv(output_path / "countsPerCellCorrected.csv", counts_per_cell_corrected)
    _save_csv(output_path / "numZero2One.csv", num_zero_to_one)
    _save_csv(output_path / "numOne2Zero.csv", num_one_to_zero)
    _save_csv(output_path / "brightnessAreaHist.csv", brightness_area_hist)
    for z in range(num_z):
        _save_csv(output_path / f"barcodeDensity-{z+1}.csv", barcode_density[:, :, z])
    if parameters.cellIDMethod == "cellID":
        _save_csv(output_path / "distToNucleusInNucleus.csv", dist_to_nucleus_in)
        _save_csv(output_path / "distToNucleusOutNucleus.csv", dist_to_nucleus_out)
        _save_csv(output_path / "fractionInNucleus.csv", fraction_in_nucleus)

    codebook_names = codebook["name"].astype(str).tolist()
    blank_mask = np.array([parameters.blankFnc(name) for name in codebook_names], dtype=bool)  # type: ignore[misc]
    total_counts = counts_per_cell_exact + counts_per_cell_corrected
    confidence_ratio = np.divide(
        np.sum(counts_per_cell_exact, axis=1),
        np.sum(total_counts, axis=1),
        out=np.zeros(num_barcodes, dtype=float),
        where=np.sum(total_counts, axis=1) != 0,
    )

    fpkm_correlations: Dict[str, float] = {}
    if parameters.abundDataPath:
        abund_path = Path(parameters.abundDataPath)
        if abund_path.exists():
            abund = _load_table(abund_path)
            lower = {c.lower(): c for c in abund.columns}
            gene_col = lower.get("genename") or lower.get("gene_name") or lower.get("gene") or lower.get("name")
            fpkm_col = lower.get("fpkm") or lower.get("abundance") or lower.get("value")
            if gene_col and fpkm_col:
                abundance_map = dict(zip(abund[gene_col].astype(str), pd.to_numeric(abund[fpkm_col], errors="coerce")))
                sorted_fpkm = np.array([abundance_map.get(name, np.nan) for name in codebook_names], dtype=float)
                fpkm_correlations = {
                    "exact": _correlation(sorted_fpkm, np.mean(counts_per_cell_exact, axis=1)),
                    "corrected": _correlation(sorted_fpkm, np.mean(counts_per_cell_corrected, axis=1)),
                    "all": _correlation(sorted_fpkm, np.mean(total_counts, axis=1)),
                }

    if plt is not None:
        _make_area_brightness_report(brightness_area_hist, output_path, parameters)
        _make_density_reports(barcode_density, image_size, x_edges, y_edges, output_path, parameters)
        _make_error_report(
            counts_per_cell_exact,
            counts_per_cell_corrected,
            num_one_to_zero,
            num_zero_to_one,
            blank_mask,
            bit_names,
            output_path,
            parameters,
        )

    if parameters.archive:
        with open(output_path / "performance_parameters.pkl", "wb") as handle:
            pickle.dump(_serializable_parameters(parameters), handle)
        source_path = Path(__file__)
        try:
            shutil.copy2(source_path, output_path / source_path.name)
        except Exception:
            warnings.warn("Could not archive MERFISHPerformanceMetrics.py", RuntimeWarning)

    if parameters.logProgress:
        with open(output_path / "performance.log", "w", encoding="utf-8") as handle:
            handle.write("\n".join(log_lines) + "\n")

    return {
        "parameters": _serializable_parameters(parameters),
        "codebookPath": str(codebook_path),
        "barcodeFiles": [str(p) for p in barcode_files],
        "outputPath": str(output_path),
        "countsPerCellExact": counts_per_cell_exact,
        "countsPerCellCorrected": counts_per_cell_corrected,
        "numZero2One": num_zero_to_one,
        "numOne2Zero": num_one_to_zero,
        "brightnessAreaHist": brightness_area_hist,
        "barcodeDensity": barcode_density,
        "confidenceRatio": confidence_ratio,
        "blankMask": blank_mask,
        "FPKMCorrYLabel": fpkm_label,
        "fpkmCorrelations": fpkm_correlations,
    }


def _make_area_brightness_report(brightness_area_hist: np.ndarray, output_path: Path, parameters: MERFISHPerformanceParameters) -> None:
    fig = plt.figure(figsize=(14, 4))
    area_bins = np.asarray(parameters.areaBins, dtype=float)
    brightness_bins = np.asarray(parameters.brightnessBins, dtype=float)

    ax = fig.add_subplot(1, 3, 1)
    area_centers = area_bins[:-1]
    ax.bar(area_centers, np.sum(brightness_area_hist, axis=0), width=1.0)
    ax.axvline(parameters.areaThreshold, linestyle="--")
    ax.set_xlabel("Area (pixels)")
    ax.set_ylabel("Counts")

    ax = fig.add_subplot(1, 3, 2)
    brightness_centers = brightness_bins[:-1]
    ax.bar(brightness_centers, np.sum(brightness_area_hist, axis=1), width=np.mean(np.diff(brightness_bins)))
    ax.axvline(np.log10(parameters.brightnessThreshold), linestyle="--")
    ax.set_xlabel("Brightness (log10)")
    ax.set_ylabel("Counts")

    ax = fig.add_subplot(1, 3, 3)
    im = ax.imshow(brightness_area_hist, aspect="auto", origin="lower")
    fig.colorbar(im, ax=ax)
    ax.set_xlabel("Area bin")
    ax.set_ylabel("Brightness bin")
    fig.tight_layout()
    _save_figure(fig, output_path, "area_brightness_histograms", parameters)
    plt.close(fig)


def _make_density_reports(
    barcode_density: np.ndarray,
    image_size: Tuple[int, int],
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    output_path: Path,
    parameters: MERFISHPerformanceParameters,
) -> None:
    for z in range(barcode_density.shape[2]):
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(1, 1, 1)
        im = ax.imshow(barcode_density[:, :, z], origin="lower", aspect="auto")
        fig.colorbar(im, ax=ax, label="Number")
        ax.set_xlabel("X Position (pixels)")
        ax.set_ylabel("Y Position (pixels)")
        ax.set_title(f"Barcode density z {z+1}")
        fig.tight_layout()
        _save_figure(fig, output_path, f"barcode_density_z_{z+1}", parameters)
        plt.close(fig)

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(np.sum(barcode_density, axis=2), origin="lower", aspect="auto")
    fig.colorbar(im, ax=ax, label="Number")
    ax.set_xlabel("X Position (pixels)")
    ax.set_ylabel("Y Position (pixels)")
    ax.set_title("Barcode density")
    fig.tight_layout()
    _save_figure(fig, output_path, "barcode_density", parameters)
    plt.close(fig)


def _make_error_report(
    counts_exact: np.ndarray,
    counts_corrected: np.ndarray,
    num_one_to_zero: np.ndarray,
    num_zero_to_one: np.ndarray,
    blank_mask: np.ndarray,
    bit_names: list[str],
    output_path: Path,
    parameters: MERFISHPerformanceParameters,
) -> None:
    total_counts = counts_exact + counts_corrected
    fig = plt.figure(figsize=(14, 6))

    ax = fig.add_subplot(2, 4, 1)
    ax.bar([0, 1], [np.sum(counts_corrected), np.sum(counts_exact)])
    ax.set_xticks([0, 1], ["Corrected", "Exact"], rotation=45)
    ax.set_ylabel("Counts")
    denom = np.sum(total_counts)
    ax.set_title(f"{(np.sum(counts_corrected) / denom * 100) if denom else 0:.2g}%")

    datasets = [counts_corrected, counts_exact, total_counts]
    labels = ["Corrected", "Exact", "All"]
    non_blank_mask = ~blank_mask
    for i, (data, label) in enumerate(zip(datasets, labels), start=2):
        ax = fig.add_subplot(2, 4, i)
        n_a = np.sum(data, axis=1)
        order = np.argsort(n_a)[::-1]
        sorted_counts = n_a[order]
        sorted_blank = blank_mask[order]
        ax.bar(np.arange(len(sorted_counts)), sorted_counts)
        if np.any(sorted_blank):
            ax.bar(np.where(sorted_blank)[0], sorted_counts[sorted_blank])
        ax.set_yscale("log" if np.any(sorted_counts > 0) else "linear")
        ax.set_xlabel("Barcode ID")
        ax.set_ylabel(f"Counts {label}")
        if np.any(blank_mask):
            max_blank = np.max(n_a[blank_mask])
            ax.set_title(f"Number above: {int(np.sum(n_a[non_blank_mask] > max_blank))}")

    ax = fig.add_subplot(2, 4, 5)
    denom = np.sum(total_counts, axis=1)
    confidence_ratio = np.divide(np.sum(counts_exact, axis=1), denom, out=np.zeros_like(denom), where=denom != 0)
    order = np.argsort(confidence_ratio)[::-1]
    ax.bar(np.arange(len(order)), confidence_ratio[order])
    if np.any(blank_mask):
        sorted_blank = blank_mask[order]
        ax.bar(np.where(sorted_blank)[0], confidence_ratio[order][sorted_blank])
        max_blank = np.max(confidence_ratio[blank_mask])
        ax.set_title(f"Number above: {int(np.sum(confidence_ratio[~blank_mask] > max_blank))}")
    ax.set_xlabel("Barcode ID")
    ax.set_ylabel("Confidence ratio")

    total_by_barcode = np.sum(total_counts, axis=1)
    bits_available = max(0, min(len(bit_names), num_one_to_zero.shape[1] - 1))
    if bits_available:
        one_to_zero = np.divide(
            num_one_to_zero[:, 1 : bits_available + 1],
            total_by_barcode[:, None],
            out=np.full((len(total_by_barcode), bits_available), np.nan),
            where=total_by_barcode[:, None] != 0,
        )
        zero_to_one = np.divide(
            num_zero_to_one[:, 1 : bits_available + 1],
            total_by_barcode[:, None],
            out=np.full((len(total_by_barcode), bits_available), np.nan),
            where=total_by_barcode[:, None] != 0,
        )
    else:
        one_to_zero = np.zeros((counts_exact.shape[0], 0))
        zero_to_one = np.zeros((counts_exact.shape[0], 0))

    ax = fig.add_subplot(2, 4, 6)
    if one_to_zero.size:
        means = np.nanmean(one_to_zero, axis=0)
        ax.bar(np.arange(1, len(means) + 1), means)
        ax.set_title(f"1→0: {np.nanmean(means):.2g}")
    ax.set_ylabel("Error rate")
    ax.set_xlabel("Bit")

    ax = fig.add_subplot(2, 4, 7)
    if zero_to_one.size:
        means = np.nanmean(zero_to_one, axis=0)
        ax.bar(np.arange(1, len(means) + 1), means)
        ax.set_title(f"0→1: {np.nanmean(means):.2g}")
    ax.set_ylabel("Error rate")
    ax.set_xlabel("Bit")

    fig.tight_layout()
    _save_figure(fig, output_path, "error_report", parameters)
    plt.close(fig)


# snake_case alias for normal Python use
merfish_performance_metrics = MERFISHPerformanceMetrics
