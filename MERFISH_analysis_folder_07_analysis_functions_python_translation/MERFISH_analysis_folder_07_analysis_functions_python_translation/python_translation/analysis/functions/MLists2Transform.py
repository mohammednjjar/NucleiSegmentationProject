"""
Python translation of analysis/functions/MLists2Transform.m from
ZhuangLab/MERFISH_analysis.

Purpose/use:
    Build one geometric transform per frame by matching points between a
    reference molecule list and a moving molecule list, optionally apply the
    transform to update moving-list xc/yc coordinates, and return residuals
    and control-point indices.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Iterable, List, Tuple
import copy
import warnings

import numpy as np

try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover - fallback path
    cKDTree = None

try:
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None

try:
    from .geometric_transforms import (
        GeometricTransform,
        identity_transform,
        fit_geometric_transform,
        transform_points_forward,
    )
except ImportError:  # script use without package import
    from geometric_transforms import (  # type: ignore
        GeometricTransform,
        identity_transform,
        fit_geometric_transform,
        transform_points_forward,
    )


@dataclass
class MLists2TransformParameters:
    transformationType: str = "nonreflectivesimilarity"
    polynomialOrder: int = 3
    controlPointMethod: str = "nearestNeighbor"
    distanceWeight: float = 0.25
    thetaWeight: float = 0.25
    histogramEdges: np.ndarray | None = None
    numNN: int = 10
    pairDistTolerance: float = 1.0
    applyTransform: bool = True
    ignoreFrames: bool = False
    transpose: bool = False
    debug: bool = False
    displayFraction: float = 0.05

    def __post_init__(self) -> None:
        if self.histogramEdges is None:
            self.histogramEdges = np.arange(-128, 129, 1, dtype=float)
        self.transformationType = str(self.transformationType)
        self.controlPointMethod = str(self.controlPointMethod)
        if self.transformationType not in {
            "nonreflectivesimilarity",
            "similarity",
            "affine",
            "projective",
            "polynomial",
        }:
            raise ValueError("unsupported transformationType")
        if self.controlPointMethod not in {"nearestNeighbor", "kNNDistanceHistogram"}:
            raise ValueError("unsupported controlPointMethod")
        if self.polynomialOrder < 0:
            raise ValueError("polynomialOrder must be non-negative")
        if self.distanceWeight < 0 or self.thetaWeight < 0:
            raise ValueError("distanceWeight and thetaWeight must be non-negative")
        if self.numNN <= 0:
            raise ValueError("numNN must be positive")
        if self.pairDistTolerance <= 0:
            raise ValueError("pairDistTolerance must be positive")


def _as_parameters(kwargs: Dict[str, Any]) -> MLists2TransformParameters:
    allowed = set(MLists2TransformParameters.__dataclass_fields__.keys())
    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(f"unknown parameter(s): {sorted(unknown)}")
    return MLists2TransformParameters(**kwargs)


def _extract_array(m_list: Any, field: str) -> np.ndarray:
    if pd is not None and isinstance(m_list, pd.DataFrame):
        if field not in m_list.columns:
            raise ValueError(f"molecule list is missing field/column {field!r}")
        return m_list[field].to_numpy()
    if isinstance(m_list, dict):
        if field not in m_list:
            raise ValueError(f"molecule list is missing field {field!r}")
        return np.asarray(m_list[field])
    if hasattr(m_list, field):
        return np.asarray(getattr(m_list, field))
    raise ValueError(f"molecule list is missing field {field!r}")


def _copy_mlist(m_list: Any) -> Any:
    if pd is not None and isinstance(m_list, pd.DataFrame):
        return m_list.copy()
    return copy.deepcopy(m_list)


def _set_array(m_list: Any, field: str, values: np.ndarray, mask: np.ndarray | None = None) -> None:
    if pd is not None and isinstance(m_list, pd.DataFrame):
        if field not in m_list.columns:
            warnings.warn(
                "Transformed points were not stored because the xc/yc fields were missing from mList",
                RuntimeWarning,
            )
            return
        if field in m_list.columns:
            m_list[field] = m_list[field].astype(float)
        if mask is None:
            m_list.loc[:, field] = np.asarray(values, dtype=float)
        else:
            m_list.loc[np.asarray(mask), field] = np.asarray(values, dtype=float)
        return
    if isinstance(m_list, dict):
        if field not in m_list:
            warnings.warn(
                "Transformed points were not stored because the xc/yc fields were missing from mList",
                RuntimeWarning,
            )
            return
        arr = np.asarray(m_list[field]).copy()
        if mask is None:
            arr[...] = values
        else:
            arr[np.asarray(mask)] = values
        m_list[field] = arr
        return
    if not hasattr(m_list, field):
        warnings.warn(
            "Transformed points were not stored because the xc/yc fields were missing from mList",
            RuntimeWarning,
        )
        return
    arr = np.asarray(getattr(m_list, field)).copy()
    if mask is None:
        arr[...] = values
    else:
        arr[np.asarray(mask)] = values
    setattr(m_list, field, arr)


def _points_for_frame(m_list: Any, frame: int | None, transpose: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(_extract_array(m_list, "x"), dtype=float)
    y = np.asarray(_extract_array(m_list, "y"), dtype=float)
    frames = np.asarray(_extract_array(m_list, "frame"))
    if frame is None:
        mask = np.ones_like(frames, dtype=bool)
    else:
        mask = frames == frame
    if transpose:
        points = np.column_stack([x[mask].reshape(-1), y[mask].reshape(-1)])
    else:
        points = np.column_stack([x[mask], y[mask]])
    return points.astype(float), mask


def _nearest_neighbors(moving: np.ndarray, reference: np.ndarray, k: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    if moving.size == 0 or reference.size == 0:
        return np.empty((0,), dtype=int), np.empty((0,), dtype=float)
    k = max(1, min(int(k), moving.shape[0]))
    if cKDTree is not None:
        distances, indices = cKDTree(moving).query(reference, k=k)
        return np.asarray(indices), np.asarray(distances, dtype=float)
    diff = reference[:, None, :] - moving[None, :, :]
    distances_all = np.sqrt(np.sum(diff * diff, axis=2))
    if k == 1:
        indices = np.argmin(distances_all, axis=1)
        return indices, distances_all[np.arange(reference.shape[0]), indices]
    indices = np.argsort(distances_all, axis=1)[:, :k]
    row = np.arange(reference.shape[0])[:, None]
    return indices, distances_all[row, indices]


def _select_control_points_nearest(
    ref_points: np.ndarray,
    mov_points: np.ndarray,
    distance_weight: float,
    theta_weight: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    idx, dist = _nearest_neighbors(mov_points, ref_points, k=1)
    idx = np.asarray(idx, dtype=int).reshape(-1)
    dist = np.asarray(dist, dtype=float).reshape(-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = np.arcsin(np.divide(mov_points[idx, 1] - ref_points[:, 1], dist))
    theta[~np.isfinite(theta)] = 0.0

    if not np.all(dist == 0):
        d_bounds = np.median(dist) + np.std(dist) * distance_weight * np.array([-1.0, 1.0])
        theta_bounds = np.median(theta) + np.std(theta) * theta_weight * np.array([-1.0, 1.0])
        keep = np.where(
            (dist >= d_bounds[0])
            & (dist <= d_bounds[1])
            & (theta >= theta_bounds[0])
            & (theta <= theta_bounds[1])
        )[0]
    else:
        keep = np.arange(len(idx), dtype=int)
    return idx, dist, keep


def _select_control_points_histogram(
    ref_points: np.ndarray,
    mov_points: np.ndarray,
    histogram_edges: np.ndarray,
    num_nn: int,
    pair_dist_tolerance: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    idx_matrix, _ = _nearest_neighbors(mov_points, ref_points, k=num_nn)
    idx_matrix = np.asarray(idx_matrix, dtype=int)
    if idx_matrix.ndim == 1:
        idx_matrix = idx_matrix[:, None]
    k_actual = idx_matrix.shape[1]
    idx1 = idx_matrix.reshape(-1)
    idx2 = np.repeat(np.arange(ref_points.shape[0], dtype=int), k_actual)[: len(idx1)]
    diff_x = mov_points[idx1, 0] - ref_points[idx2, 0]
    diff_y = mov_points[idx1, 1] - ref_points[idx2, 1]
    hist, x_edges, y_edges = np.histogram2d(diff_x, diff_y, bins=[histogram_edges, histogram_edges])
    max_ind = np.unravel_index(np.argmax(hist), hist.shape)
    max_value = hist[max_ind]
    offset = np.array([histogram_edges[max_ind[0]], histogram_edges[max_ind[1]]], dtype=float)

    shifted_moving = mov_points - offset
    idx, dist = _nearest_neighbors(shifted_moving, ref_points, k=1)
    idx = np.asarray(idx, dtype=int).reshape(-1)
    dist = np.asarray(dist, dtype=float).reshape(-1)
    step = float(np.mean(np.diff(histogram_edges))) if len(histogram_edges) > 1 else 1.0
    keep = np.where(dist <= pair_dist_tolerance * step)[0]

    if max_value < np.mean(hist) + 2 * np.std(hist):
        warnings.warn(
            "The maximum kNN offset is <2 standard deviations from the mean; control point identification may be inaccurate",
            RuntimeWarning,
        )
    return idx, dist, keep


def MLists2Transform(refList: Any, mList: Any, **kwargs: Any):
    """
    Translate MATLAB MLists2Transform.

    Parameters
    ----------
    refList, mList:
        dict-like objects, pandas DataFrames, or objects with x, y, frame
        arrays. mList may contain xc and yc arrays; these are updated when
        applyTransform=True.

    Keyword arguments match the MATLAB parameter names.

    Returns
    -------
    tforms, updated_mList, residuals, inds, parameters
    """
    for field in ("x", "y", "frame"):
        _extract_array(refList, field)
        _extract_array(mList, field)
    parameters = _as_parameters(kwargs)
    moving_out = _copy_mlist(mList)

    ref_frames = np.asarray(_extract_array(refList, "frame"))
    if parameters.ignoreFrames:
        frames: List[int | None] = [None]
    else:
        frames = list(range(1, int(np.max(ref_frames)) + 1)) if ref_frames.size else []

    tforms: List[GeometricTransform] = [identity_transform() for _ in frames]
    residuals: List[np.ndarray] = []
    inds: List[np.ndarray] = []

    for frame_index, frame in enumerate(frames):
        ref_points, _ = _points_for_frame(refList, frame, transpose=parameters.transpose)
        mov_points, mov_mask = _points_for_frame(mList, frame, transpose=parameters.transpose)

        if ref_points.size == 0:
            warnings.warn("No molecules in reference list. Using default transformation.", RuntimeWarning)
            residuals.append(np.zeros((0, 4), dtype=float))
            inds.append(np.zeros((0, 2), dtype=int))
            continue
        if mov_points.size == 0:
            warnings.warn("No molecules in moving list. Using default transformation.", RuntimeWarning)
            residuals.append(np.zeros((0, 4), dtype=float))
            inds.append(np.zeros((0, 2), dtype=int))
            continue

        if parameters.controlPointMethod == "nearestNeighbor":
            idx, distances, points_to_keep = _select_control_points_nearest(
                ref_points,
                mov_points,
                parameters.distanceWeight,
                parameters.thetaWeight,
            )
        else:
            idx, distances, points_to_keep = _select_control_points_histogram(
                ref_points,
                mov_points,
                np.asarray(parameters.histogramEdges, dtype=float),
                parameters.numNN,
                parameters.pairDistTolerance,
            )

        try:
            selected_moving = mov_points[idx[points_to_keep], :]
            selected_fixed = ref_points[points_to_keep, :]
            tforms[frame_index] = fit_geometric_transform(
                selected_moving,
                selected_fixed,
                transformation_type=parameters.transformationType,
                polynomial_order=parameters.polynomialOrder,
            )
        except Exception as exc:
            warnings.warn(f"Did not find sufficient control points. Using default transformation. Reason: {exc}", RuntimeWarning)
            tforms[frame_index] = identity_transform()

        inds.append(np.column_stack([points_to_keep, idx[points_to_keep]]).astype(int))

        if parameters.applyTransform:
            moved_points = transform_points_forward(tforms[frame_index], mov_points)
            _set_array(moving_out, "xc", moved_points[:, 0], mask=None if parameters.ignoreFrames else mov_mask)
            _set_array(moving_out, "yc", moved_points[:, 1], mask=None if parameters.ignoreFrames else mov_mask)
            if points_to_keep.size:
                local_residuals = np.column_stack([
                    moved_points[idx[points_to_keep], :] - ref_points[points_to_keep, :],
                    mov_points[idx[points_to_keep], :],
                ])
            else:
                local_residuals = np.zeros((0, 4), dtype=float)
            residuals.append(local_residuals)
        else:
            residuals.append(np.zeros((0, 4), dtype=float))

    if parameters.ignoreFrames:
        tforms_out = tforms[0] if tforms else identity_transform()
        residuals_out = residuals[0] if residuals else np.zeros((0, 4), dtype=float)
        inds_out = inds[0] if inds else np.zeros((0, 2), dtype=int)
    else:
        tforms_out = tforms
        residuals_out = residuals
        inds_out = inds

    return tforms_out, moving_out, residuals_out, inds_out, asdict(parameters)


# snake_case alias for normal Python use
mlists_to_transform = MLists2Transform
