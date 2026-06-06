
"""Shared utilities for the translated `analysis/deprecated` MATLAB folder.

This module replaces MATLAB struct-array access, containers.Map key conversion,
nearest-neighbor search, affine/rigid transforms, and simple word/molecule-list
helpers used by the old MERFISH analysis scripts.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import os
import pickle
import re
import uuid
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

Record = dict[str, Any]


def as_record(obj: Any) -> Record:
    if isinstance(obj, Mapping):
        return dict(obj)
    return dict(vars(obj))


def ensure_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return list(value)
        if value.ndim == 0:
            return [value.item()]
        return list(value)
    return [value]


def get_field(obj: Any, name: str, default: Any = None) -> Any:
    if isinstance(obj, Mapping):
        return obj.get(name, default)
    return getattr(obj, name, default)


def set_field(obj: Any, name: str, value: Any) -> None:
    if isinstance(obj, dict):
        obj[name] = value
    else:
        setattr(obj, name, value)


def has_field(obj: Any, name: str) -> bool:
    if isinstance(obj, Mapping):
        return name in obj
    return hasattr(obj, name)


def fields(obj: Any) -> list[str]:
    if isinstance(obj, Mapping):
        return list(obj.keys())
    return list(vars(obj).keys())


def remove_fields(obj: Any, names: Iterable[str]) -> Any:
    for name in list(names):
        if isinstance(obj, dict):
            obj.pop(name, None)
        elif hasattr(obj, name):
            delattr(obj, name)
    return obj


def copy_selected_fields(src: Any, dst: Any, src_names: Sequence[str], dst_names: Sequence[str]) -> None:
    for src_name, dst_name in zip(src_names, dst_names):
        set_field(dst, dst_name, get_field(src, src_name))


def parse_parameters(defaults: Mapping[str, Any] | None = None, parameters: Mapping[str, Any] | None = None, **kwargs: Any) -> Record:
    out: Record = dict(defaults or {})
    if parameters is not None:
        out.update(dict(parameters))
    out.update(kwargs)
    return out


def uuid_string() -> str:
    return str(uuid.uuid4())


def bits_to_binstr(bits: Sequence[int] | np.ndarray | str) -> str:
    if isinstance(bits, str):
        return re.sub(r"\s+", "", bits)
    return "".join(str(int(x)) for x in np.asarray(bits).astype(int).ravel().tolist())


def bits_to_int(bits: Sequence[int] | np.ndarray | str) -> int:
    clean = bits_to_binstr(bits)
    return int(clean, 2) if clean else 0


def key_converter(key_type: str) -> Callable[[Sequence[int] | np.ndarray | str], Any]:
    key_type = str(key_type)
    if key_type == "int":
        return bits_to_int
    if key_type == "binStr":
        return bits_to_binstr
    raise ValueError("key_type must be 'int' or 'binStr'")


def get_array_field(obj: Any, name: str, default: Any = None, dtype: Any = float) -> np.ndarray:
    value = get_field(obj, name, default)
    if value is None:
        return np.asarray([], dtype=dtype)
    return np.asarray(value, dtype=dtype).ravel()


def mlist_xy(mlist: Any, x_field: str = "xc", y_field: str = "yc") -> np.ndarray:
    x = get_array_field(mlist, x_field)
    y = get_array_field(mlist, y_field)
    if len(x) != len(y):
        n = min(len(x), len(y))
        x, y = x[:n], y[:n]
    if len(x) == 0:
        return np.zeros((0, 2), dtype=float)
    return np.column_stack([x, y]).astype(float)


def filter_mlist(mlist: Any, frame: Sequence[int] | None = None) -> Record:
    rec = as_record(mlist)
    if frame is None or "frame" not in rec:
        return rec
    allowed = set(int(x) for x in frame)
    fr = np.asarray(rec.get("frame"), dtype=int).ravel()
    mask = np.array([int(x) in allowed for x in fr], dtype=bool)
    out: Record = {}
    for key, val in rec.items():
        arr = np.asarray(val)
        if arr.ndim >= 1 and arr.shape[0] == mask.shape[0]:
            out[key] = arr[mask].tolist()
        else:
            out[key] = val
    return out


def nearest_neighbors(query: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    query = np.asarray(query, dtype=float).reshape(-1, 2)
    reference = np.asarray(reference, dtype=float).reshape(-1, 2)
    if len(query) == 0 or len(reference) == 0:
        return np.full(len(query), -1, dtype=int), np.full(len(query), np.inf)
    diff = query[:, None, :] - reference[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    idx = np.argmin(dist, axis=1)
    return idx.astype(int), dist[np.arange(len(query)), idx]


def mutual_nearest_pairs(reference: np.ndarray, moving: np.ndarray, max_distance: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference = np.asarray(reference, dtype=float).reshape(-1, 2)
    moving = np.asarray(moving, dtype=float).reshape(-1, 2)
    if len(reference) == 0 or len(moving) == 0:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    idx_ref, dist_ref = nearest_neighbors(moving, reference)
    idx_mov, _dist_mov = nearest_neighbors(reference, moving)
    moving_ids: list[int] = []
    reference_ids: list[int] = []
    distances: list[float] = []
    for moving_i, ref_i in enumerate(idx_ref):
        if ref_i >= 0 and idx_mov[ref_i] == moving_i and dist_ref[moving_i] <= max_distance:
            reference_ids.append(int(ref_i))
            moving_ids.append(int(moving_i))
            distances.append(float(dist_ref[moving_i]))
    return np.asarray(reference_ids, dtype=int), np.asarray(moving_ids, dtype=int), np.asarray(distances, dtype=float)


@dataclass
class AffineTransform:
    """Affine transform with a 3x3 homogeneous matrix mapping reference to moving.

    MATLAB code usually calls `tforminv(tform, moving_x, moving_y)` to put moving
    image coordinates into the reference frame. Here `inverse_points` implements
    that same operation.
    """
    matrix: np.ndarray

    @staticmethod
    def identity() -> "AffineTransform":
        return AffineTransform(np.eye(3, dtype=float))

    def forward_points(self, points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=float).reshape(-1, 2)
        if len(pts) == 0:
            return np.zeros((0, 2), dtype=float)
        hom = np.column_stack([pts, np.ones(len(pts))])
        out = hom @ self.matrix.T
        return out[:, :2]

    def inverse_points(self, points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=float).reshape(-1, 2)
        if len(pts) == 0:
            return np.zeros((0, 2), dtype=float)
        inv = np.linalg.inv(self.matrix)
        hom = np.column_stack([pts, np.ones(len(pts))])
        out = hom @ inv.T
        return out[:, :2]

    def compose_after(self, previous: "AffineTransform") -> "AffineTransform":
        return AffineTransform(self.matrix @ previous.matrix)


def estimate_rigid_transform(reference: np.ndarray, moving: np.ndarray) -> AffineTransform:
    """Estimate a rotation+translation transform mapping reference points to moving points."""
    ref = np.asarray(reference, dtype=float).reshape(-1, 2)
    mov = np.asarray(moving, dtype=float).reshape(-1, 2)
    if ref.shape != mov.shape or len(ref) == 0:
        raise ValueError("reference and moving must have the same non-zero shape")
    if len(ref) == 1:
        shift = mov[0] - ref[0]
        mat = np.eye(3, dtype=float)
        mat[:2, 2] = shift
        return AffineTransform(mat)
    ref_centroid = ref.mean(axis=0)
    mov_centroid = mov.mean(axis=0)
    ref0 = ref - ref_centroid
    mov0 = mov - mov_centroid
    cov = ref0.T @ mov0
    u, _s, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    trans = mov_centroid - ref_centroid @ rot.T
    mat = np.eye(3, dtype=float)
    mat[:2, :2] = rot
    mat[:2, 2] = trans
    return AffineTransform(mat)


def compose_transforms(transforms: Sequence[AffineTransform | None]) -> AffineTransform:
    current = AffineTransform.identity()
    for tform in transforms:
        if tform is not None:
            current = AffineTransform(tform.matrix @ current.matrix)
    return current


def create_words_structure(num_words: int, num_hybs: int) -> list[Record]:
    out: list[Record] = []
    for _ in range(int(num_words)):
        out.append({
            "measuredCodeword": np.zeros(int(num_hybs), dtype=bool),
            "codeword": np.zeros(int(num_hybs), dtype=bool),
            "mListInds": [],
            "wordCentroidX": np.nan,
            "wordCentroidY": np.nan,
            "intCodeword": 0,
            "numOnBits": 0,
            "measuredOnBits": [],
            "onBits": [],
            "geneName": "",
            "isExactMatch": False,
            "isCorrectedMatch": False,
        })
    return out


def stable_unique_rows(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float).reshape(-1, 2)
    seen: set[tuple[Any, Any]] = set()
    rows: list[np.ndarray] = []
    idxs: list[int] = []
    for idx, row in enumerate(values):
        key = tuple(None if np.isnan(x) else round(float(x), 12) for x in row)
        if key in seen:
            continue
        seen.add(key)
        rows.append(row.copy())
        idxs.append(idx)
    if rows:
        return np.vstack(rows), np.asarray(idxs, dtype=int)
    return np.zeros((0, 2), dtype=float), np.array([], dtype=int)


def read_fasta(path: str | os.PathLike[str]) -> list[Record]:
    records: list[Record] = []
    header = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"Header": header, "Sequence": "".join(chunks)})
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        records.append({"Header": header, "Sequence": "".join(chunks)})
    return records


def simple_codebook_to_map(codebook: Sequence[Mapping[str, Any]], key_type: str = "binStr") -> dict[Any, Any]:
    conv = key_converter(key_type)
    out: dict[Any, Any] = {}
    for rec in codebook:
        seq = rec.get("Sequence", rec.get("sequence", rec.get("codeword", "")))
        name = rec.get("Header", rec.get("name", rec.get("geneName", "")))
        out[conv(seq)] = name
    return out


def parse_value(text: str) -> Any:
    text = str(text).strip()
    if text == "":
        return text
    for caster in (int, float):
        try:
            return caster(text)
        except ValueError:
            continue
    low = text.lower()
    if low in {"true", "false"}:
        return low == "true"
    return text


def read_master_molecule_list(path: str | os.PathLike[str], compact: bool = True, transpose: bool = True, verbose: bool = False) -> Record:
    """Read a molecule list from common Python-friendly files.

    The old MATLAB script reads STORM `.bin` files through matlab-storm. This
    translation supports CSV/TSV/JSON/NPY/NPZ/Pickle molecule lists directly.
    """
    p = Path(path)
    ext = p.suffix.lower()
    if ext in {".csv", ".txt", ".tsv"}:
        delimiter = "\t" if ext == ".tsv" else ","
        with p.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            columns: dict[str, list[Any]] = {}
            for row in reader:
                for key, value in row.items():
                    columns.setdefault(key, []).append(parse_value(value))
        return columns
    if ext == ".json":
        data = json.loads(p.read_text(encoding="utf-8"))
        if isinstance(data, list):
            keys = sorted({k for row in data if isinstance(row, Mapping) for k in row.keys()})
            return {key: [row.get(key) for row in data] for key in keys}
        if isinstance(data, Mapping):
            return dict(data)
    if ext == ".npy":
        arr = np.load(p, allow_pickle=True)
        if arr.dtype.names:
            return {name: arr[name].tolist() for name in arr.dtype.names}
        return {"array": arr.tolist()}
    if ext == ".npz":
        data = np.load(p, allow_pickle=True)
        return {key: data[key].tolist() for key in data.files}
    if ext in {".pkl", ".pickle"}:
        with p.open("rb") as handle:
            obj = pickle.load(handle)
        return as_record(obj) if not isinstance(obj, list) else {"records": obj}
    raise ValueError(f"Unsupported molecule-list extension for Python reader: {p.suffix}")


def build_image_data_structures(data_path: str | os.PathLike[str], parameters: Mapping[str, Any] | None = None) -> list[Record]:
    """Directory scanner corresponding to MATLAB BuildImageDataStructures.

    It extracts `hybNum`, `cellNum`, `isFiducial`, `movieType`, and `binType`
    from file names using flexible keyword rules so AnalyzeMERFISH can run on
    Python-readable molecule lists.
    """
    params = dict(parameters or {})
    file_ext = str(params.get("fileExt", "bin")).lstrip(".")
    image_tag = str(params.get("imageTag", "STORM"))
    root = Path(data_path)
    files = sorted(p for p in root.rglob(f"*.{file_ext}") if p.is_file())
    records: list[Record] = []
    for path in files:
        stem = path.stem
        lower = stem.lower()
        hyb_match = re.search(r"(?:hyb|h)(\d+)", lower)
        cell_match = re.search(r"(?:cell|fov|c)(\d+)", lower)
        hyb_num = int(hyb_match.group(1)) if hyb_match else 1
        cell_num = int(cell_match.group(1)) if cell_match else 1
        is_fid = "fid" in lower or "bead" in lower or "c2" in lower
        bin_type = str(params.get("fiducialMListType" if is_fid else "imageMListType", "list" if is_fid else "alist"))
        records.append({
            "name": stem,
            "filePath": str(path),
            "movieType": image_tag,
            "hybNum": hyb_num,
            "cellNum": cell_num,
            "isFiducial": bool(is_fid),
            "binType": bin_type,
        })
    return records


def transfer_info_file_fields(records: Sequence[Any], parameters: Mapping[str, Any] | None = None) -> list[Any]:
    """Python equivalent hook for TransferInfoFileFields.

    It preserves input records and fills missing stage/image metadata with stable defaults.
    """
    out = ensure_list(records)
    for rec in out:
        for key, value in {"Stage_X": 0.0, "Stage_Y": 0.0, "focusLockQuality": np.nan}.items():
            if not has_field(rec, key):
                set_field(rec, key, value)
    return out


def get_report_visibility(reports_to_generate: Any, name: str) -> bool:
    for item in ensure_list(reports_to_generate):
        if isinstance(item, (list, tuple)) and len(item) >= 2 and item[0] == name:
            return bool(item[1])
        if item == name:
            return True
    return False
