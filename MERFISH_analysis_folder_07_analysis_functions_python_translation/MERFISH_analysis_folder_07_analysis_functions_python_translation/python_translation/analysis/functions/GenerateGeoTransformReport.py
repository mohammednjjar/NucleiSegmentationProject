"""
Python translation of analysis/functions/GenerateGeoTransformReport.m.

Purpose/use:
    Summarize geometric-transform residual errors, spatial dependence of
    residuals, and transform offset/scale/angle. Optionally create matplotlib
    figures matching the MATLAB report categories.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Iterable, List, Sequence, Tuple
import os
import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

try:
    from .geometric_transforms import transform_points_forward
except ImportError:
    from geometric_transforms import transform_points_forward  # type: ignore


@dataclass
class GeoTransformReportParameters:
    reportsToGenerate: Tuple[Tuple[str, str, Tuple[int, int, int, int]], ...] = (
        ("residualTransformError", "on", (1, 1, 1700, 400)),
        ("residualTransformErrorByPosition", "on", (1, 1, 1700, 400)),
        ("transformSummary", "on", (1, 1, 1230, 400)),
    )
    edges: Tuple[np.ndarray, np.ndarray] | None = None
    saveAndClose: bool = False
    overwrite: bool = True
    formats: Tuple[str, ...] = ("png", "fig")
    savePath: str | None = None
    makeFigures: bool = True

    def __post_init__(self) -> None:
        if self.edges is None:
            self.edges = (np.arange(1, 513, 25, dtype=float), np.arange(1, 513, 25, dtype=float))
        else:
            self.edges = (np.asarray(self.edges[0], dtype=float), np.asarray(self.edges[1], dtype=float))


def _to_object_array(items: Any) -> np.ndarray:
    if isinstance(items, np.ndarray) and items.dtype == object:
        return items
    if isinstance(items, np.ndarray):
        arr = np.empty(items.shape, dtype=object)
        for idx in np.ndindex(items.shape):
            arr[idx] = items[idx]
        return arr
    if isinstance(items, (list, tuple)):
        if items and isinstance(items[0], np.ndarray):
            arr = np.empty((len(items),), dtype=object)
            for i, value in enumerate(items):
                arr[i] = value
            return arr
        if items and isinstance(items[0], (list, tuple)) and not hasattr(items[0], "transform_points"):
            rows = len(items)
            cols = max(len(row) for row in items)
            arr = np.empty((rows, cols), dtype=object)
            arr[:] = None
            for i, row in enumerate(items):
                for j, value in enumerate(row):
                    arr[i, j] = value
            return arr
        arr = np.empty((len(items),), dtype=object)
        for i, value in enumerate(items):
            arr[i] = value
        return arr
    arr = np.empty((1,), dtype=object)
    arr[0] = items
    return arr


def _flatten_residuals(residuals: Any) -> List[np.ndarray]:
    arr = _to_object_array(residuals)
    out: List[np.ndarray] = []
    for item in arr.flat:
        if item is None:
            out.append(np.zeros((0, 4), dtype=float))
        else:
            local = np.asarray(item, dtype=float)
            if local.size == 0:
                out.append(np.zeros((0, 4), dtype=float))
            else:
                if local.ndim != 2 or local.shape[1] < 4:
                    raise ValueError("each residual array must be N x 4 or wider")
                out.append(local[:, :4])
    return out


def _has_report(parameters: GeoTransformReportParameters, name: str) -> bool:
    return any(row[0] == name for row in parameters.reportsToGenerate)


def _safe_quantile(data: np.ndarray, q: Sequence[float]) -> np.ndarray:
    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.array([0.0, 1.0])
    return np.quantile(finite, q)


def _save_figure(fig: Any, title: str, parameters: GeoTransformReportParameters) -> None:
    if not parameters.savePath:
        return
    os.makedirs(parameters.savePath, exist_ok=True)
    base = os.path.join(parameters.savePath, title.replace(" ", "_").replace("/", "_"))
    for fmt in parameters.formats:
        if fmt.lower() == "fig":
            continue
        path = f"{base}.{fmt}"
        if os.path.exists(path) and not parameters.overwrite:
            raise FileExistsError(path)
        fig.savefig(path, bbox_inches="tight")


def GenerateGeoTransformReport(tforms: Any, residuals: Any, **kwargs: Any):
    """
    Translate MATLAB GenerateGeoTransformReport.

    Returns
    -------
    report, fig_handles, parameters
    """
    allowed = set(GeoTransformReportParameters.__dataclass_fields__.keys())
    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(f"unknown parameter(s): {sorted(unknown)}")
    parameters = GeoTransformReportParameters(**kwargs)

    tform_arr = _to_object_array(tforms)
    residual_list = _flatten_residuals(residuals)
    if len(residual_list) != tform_arr.size:
        if len(residual_list) == 1 and tform_arr.size > 1:
            residual_list = residual_list * tform_arr.size
        else:
            raise ValueError("tforms and residuals must contain the same number of elements")

    display_method = "1D" if 1 in tform_arr.shape else "2D"
    shape = tform_arr.shape
    report: dict[str, np.ndarray] = {}
    fig_handles: list[Any] = []

    if _has_report(parameters, "residualTransformError"):
        mu_x = np.full(shape, np.nan)
        mu_y = np.full(shape, np.nan)
        std_x = np.full(shape, np.nan)
        std_y = np.full(shape, np.nan)
        num_cp = np.full(shape, np.nan)
        for linear_index, local in enumerate(residual_list):
            multi_index = np.unravel_index(linear_index, shape)
            if local.size:
                mu_x[multi_index] = np.mean(local[:, 0])
                mu_y[multi_index] = np.mean(local[:, 1])
                std_x[multi_index] = np.std(local[:, 0], ddof=1) if local.shape[0] > 1 else 0.0
                std_y[multi_index] = np.std(local[:, 1], ddof=1) if local.shape[0] > 1 else 0.0
                num_cp[multi_index] = local.shape[0]
        report.update(muX=mu_x, muY=mu_y, stdX=std_x, stdY=std_y, numCP=num_cp)

    if _has_report(parameters, "residualTransformErrorByPosition"):
        non_empty = [r for r in residual_list if r.size]
        all_residuals = np.vstack(non_empty) if non_empty else np.zeros((0, 4), dtype=float)
        edges_x, edges_y = parameters.edges  # type: ignore[misc]
        err_x = np.full((len(edges_x) - 1, len(edges_y) - 1), np.nan)
        err_y = np.full_like(err_x, np.nan)
        num_values = np.zeros_like(err_x)
        if all_residuals.size:
            for i in range(len(edges_x) - 1):
                for j in range(len(edges_y) - 1):
                    mask = (
                        (all_residuals[:, 2] >= edges_x[i])
                        & (all_residuals[:, 2] < edges_x[i + 1])
                        & (all_residuals[:, 3] >= edges_y[j])
                        & (all_residuals[:, 3] < edges_y[j + 1])
                    )
                    if np.any(mask):
                        err_x[i, j] = np.mean(all_residuals[mask, 0])
                        err_y[i, j] = np.mean(all_residuals[mask, 1])
                    num_values[i, j] = int(np.sum(mask))
        report.update(errX=err_x, errY=err_y, numValues=num_values)

    if _has_report(parameters, "transformSummary"):
        offset_x = np.full(shape, np.nan)
        offset_y = np.full(shape, np.nan)
        angle = np.full(shape, np.nan)
        scale = np.full(shape, np.nan)
        unit_vector = np.array([[0.0, 0.0], [1.0, 0.0]])
        for linear_index, transform in enumerate(tform_arr.flat):
            multi_index = np.unravel_index(linear_index, shape)
            if transform is None:
                continue
            moved = transform_points_forward(transform, unit_vector)
            offset_x[multi_index] = moved[0, 0]
            offset_y[multi_index] = moved[0, 1]
            diff_vector = moved[1, :] - moved[0, :]
            local_scale = float(np.sqrt(np.sum(diff_vector * diff_vector)))
            scale[multi_index] = local_scale
            if local_scale > 0:
                cosine = float(np.clip(np.dot(diff_vector, unit_vector[1]) / local_scale, -1.0, 1.0))
                angle[multi_index] = np.degrees(np.arccos(cosine))
        report.update(offsetX=offset_x, offsetY=offset_y, scale=scale, angle=angle)

    if parameters.makeFigures and plt is not None:
        if _has_report(parameters, "residualTransformError") and {"muX", "muY", "numCP", "stdX", "stdY"}.issubset(report):
            fig = plt.figure(figsize=(14, 4))
            data = [report["muX"], report["muY"], report["numCP"], report["stdX"], report["stdY"]]
            titles = ["Residual X", "Residual Y", "Number of Points", "STD X", "STD Y"]
            q_ranges = [(0.1, 0.9), (0.1, 0.9), (0, 1), (0, 0.9), (0, 0.9)]
            for i, (local, title, q_range) in enumerate(zip(data, titles, q_ranges), start=1):
                ax = fig.add_subplot(2, 3, i)
                if display_method == "1D":
                    ax.plot(np.ravel(local))
                    ax.set_xlabel("FOV")
                    ax.set_ylabel(title)
                    limits = _safe_quantile(local, q_range)
                    if limits[1] > limits[0]:
                        ax.set_ylim(limits)
                else:
                    im = ax.imshow(local, aspect="auto")
                    fig.colorbar(im, ax=ax)
                    ax.set_xlabel("FOV")
                    ax.set_ylabel("Imaging Rounds")
                    ax.set_title(title)
            fig.tight_layout()
            _save_figure(fig, "Geometric Transformation Error Report", parameters)
            fig_handles.append(fig)
            if parameters.saveAndClose:
                plt.close(fig)

        if _has_report(parameters, "residualTransformErrorByPosition") and {"errX", "errY", "numValues"}.issubset(report):
            fig = plt.figure(figsize=(14, 4))
            data = [report["errX"], report["errY"], report["numValues"]]
            titles = ["X", "Y", "Number"]
            for i, (local, title) in enumerate(zip(data, titles), start=1):
                ax = fig.add_subplot(1, 3, i)
                im = ax.imshow(local.T, origin="lower", aspect="auto")
                fig.colorbar(im, ax=ax)
                ax.set_xlabel("Position X (Pixels)")
                ax.set_ylabel("Position Y (Pixels)")
                ax.set_title(title)
            fig.tight_layout()
            _save_figure(fig, "Geometric Transformation Error Position Dependence Report", parameters)
            fig_handles.append(fig)
            if parameters.saveAndClose:
                plt.close(fig)

        if _has_report(parameters, "transformSummary") and {"offsetX", "offsetY", "scale", "angle"}.issubset(report):
            fig = plt.figure(figsize=(10, 5))
            data = [report["offsetX"], report["offsetY"], report["scale"], report["angle"]]
            titles = ["Offset X", "Offset Y", "Scale", "Angle"]
            ylabels = ["Pixels", "Pixels", "Fraction", "Degrees"]
            for i, (local, title, ylabel) in enumerate(zip(data, titles, ylabels), start=1):
                ax = fig.add_subplot(2, 2, i)
                if display_method == "1D":
                    ax.plot(np.ravel(local))
                    ax.set_xlabel("FOV")
                    ax.set_ylabel(ylabel)
                    ax.set_title(title)
                else:
                    im = ax.imshow(local, aspect="auto")
                    fig.colorbar(im, ax=ax)
                    ax.set_xlabel("FOV")
                    ax.set_ylabel("Imaging Rounds")
                    ax.set_title(title)
            fig.tight_layout()
            _save_figure(fig, "Geometric Transformation Summary", parameters)
            fig_handles.append(fig)
            if parameters.saveAndClose:
                plt.close(fig)

    return report, fig_handles, asdict(parameters)


# snake_case alias for normal Python use
generate_geo_transform_report = GenerateGeoTransformReport
