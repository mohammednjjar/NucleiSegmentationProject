"""Python translation of analysis/classes/FoundFeature.m.

The class stores segmented MERFISH features/cells, converts pixel boundaries
into absolute coordinates, computes morphology, tests containment/overlap,
joins broken boundaries, generates masks, and exports feature tables.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple
import copy as _copy
import math
import uuid

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure, draw, morphology

try:
    from shapely.geometry import Point, Polygon, MultiPolygon
except Exception:  # Shapely is optional at import time.
    Point = None
    Polygon = None
    MultiPolygon = None


def _as_array(x: Any, dtype=float) -> np.ndarray:
    if x is None:
        return np.array([], dtype=dtype)
    arr = np.asarray(x, dtype=dtype)
    return arr


def _polygon_area(points: np.ndarray) -> float:
    points = np.asarray(points, dtype=float)
    if len(points) < 3:
        return 0.0
    x = points[:, 0]
    y = points[:, 1]
    return float(abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def _polygon_centroid(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if len(points) == 0:
        return np.array([np.nan, np.nan], dtype=float)
    if len(points) < 3:
        return np.nanmean(points, axis=0)
    x = points[:, 0]
    y = points[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    area2 = np.sum(cross)
    if abs(area2) < 1e-12:
        return np.nanmean(points, axis=0)
    cx = np.sum((x + np.roll(x, -1)) * cross) / (3.0 * area2)
    cy = np.sum((y + np.roll(y, -1)) * cross) / (3.0 * area2)
    return np.array([cx, cy], dtype=float)


def _close_boundary(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if len(points) == 0:
        return points.reshape(0, 2)
    if not np.allclose(points[0], points[-1], equal_nan=True):
        points = np.vstack([points, points[0]])
    return points


def _poly(points: np.ndarray):
    if Polygon is None or points is None or len(points) < 3:
        return None
    p = Polygon(points)
    if not p.is_valid:
        p = p.buffer(0)
    if p.is_empty:
        return None
    return p


def _interp_line(a: np.ndarray, b: np.ndarray, step: float) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = float(np.linalg.norm(a - b))
    n = max(int(round(d / max(step, 1e-12))), 0)
    if n <= 0:
        return np.empty((0, 2), dtype=float)
    xs = np.linspace(a[0], b[0], n + 2)[1:-1]
    ys = np.linspace(a[1], b[1], n + 2)[1:-1]
    return np.column_stack([xs, ys])


@dataclass
class FoundFeature:
    """Container for one segmented MERFISH feature/cell."""

    label_mat: Optional[np.ndarray] = None
    fov_id: Optional[Any] = None
    fov_center_pos: Optional[Sequence[float]] = None
    pixel_size: float = 1.0
    stage_orientation: Sequence[float] = (1.0, 1.0)
    bounding_box: Optional[Sequence[float]] = None
    z_pos: Optional[Sequence[float]] = None
    feature_label: int = 1
    verbose: bool = False
    name: str = ""
    type: str = ""

    version: str = "1.0"
    image_size: Tuple[int, int] = (0, 0)
    uID: str = field(default_factory=lambda: str(uuid.uuid4()))
    joinedUIDs: List[str] = field(default_factory=list)
    fovID: List[Any] = field(default_factory=list)
    feature_label_value: Any = None
    num_zPos: int = 0
    boundaries: List[List[np.ndarray]] = field(default_factory=list)
    abs_zPos: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    abs_boundaries: List[np.ndarray] = field(default_factory=list)
    volume: float = 0.0
    abs_volume: float = 0.0
    boundary_area: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    abs_boundary_area: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    is_broken: bool = False
    num_joined_features: int = 0
    feature_id: int = -1
    children: List[Any] = field(default_factory=list)
    parents: List[Any] = field(default_factory=list)
    metaData: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.label_mat is None:
            self.feature_label_value = self.feature_label
            if self.fov_id is not None:
                self.fovID = [self.fov_id]
            return

        labels = np.asarray(self.label_mat)
        if labels.ndim == 2:
            labels = labels[:, :, None]
        if labels.ndim != 3:
            raise ValueError("label_mat must be 2D or 3D")

        self.feature_label_value = self.feature_label
        self.image_size = tuple(int(v) for v in labels.shape[:2])
        self.fovID = [self.fov_id]
        self.num_zPos = int(labels.shape[2])
        self.abs_zPos = np.asarray(
            self.z_pos if self.z_pos is not None else np.arange(self.num_zPos), dtype=float
        )
        if self.abs_zPos.size != self.num_zPos:
            raise ValueError("z_pos length must match number of z slices")

        center = np.asarray(self.fov_center_pos if self.fov_center_pos is not None else (0.0, 0.0), dtype=float)
        orientation = np.asarray(self.stage_orientation, dtype=float)
        if orientation.size < 2:
            orientation = np.array([1.0, 1.0], dtype=float)
        if self.bounding_box is None:
            bbox = np.array([0.0, 0.0, self.image_size[1], self.image_size[0]], dtype=float)
        else:
            bbox = np.asarray(self.bounding_box, dtype=float)
        self.boundaries = [[]]
        self.abs_boundaries = []

        for z in range(self.num_zPos):
            mask = labels[:, :, z] == self.feature_label
            contours = measure.find_contours(mask.astype(float), 0.5)
            if len(contours) == 0:
                boundary_px = np.empty((0, 2), dtype=float)
            else:
                contour = max(contours, key=len)
                boundary_px = np.column_stack([contour[:, 1], contour[:, 0]])
            self.boundaries[0].append(boundary_px)
            abs_boundary = self._pixel_boundary_to_abs(boundary_px, center, orientation, bbox)
            self.abs_boundaries.append(abs_boundary)
            if len(boundary_px):
                h, w = mask.shape
                edge = (
                    np.any(boundary_px[:, 0] <= 0.5)
                    or np.any(boundary_px[:, 0] >= w - 1.5)
                    or np.any(boundary_px[:, 1] <= 0.5)
                    or np.any(boundary_px[:, 1] >= h - 1.5)
                )
                self.is_broken = bool(self.is_broken or edge)
        self.CalculateProperties()

    def _pixel_boundary_to_abs(
        self, boundary: np.ndarray, center: np.ndarray, orientation: np.ndarray, bbox: np.ndarray
    ) -> np.ndarray:
        boundary = np.asarray(boundary, dtype=float)
        if boundary.size == 0:
            return boundary.reshape(0, 2)
        pixel_um = float(self.pixel_size) / 1000.0
        x = (boundary[:, 0] - self.image_size[1] / 2.0) * pixel_um * orientation[0] + center[0]
        y = (boundary[:, 1] - self.image_size[0] / 2.0) * pixel_um * orientation[1] + center[1]
        if bbox.size >= 2:
            x = x + 0.0 * bbox[0]
            y = y + 0.0 * bbox[1]
        return np.column_stack([x, y])

    def CalculateCentroid(self) -> np.ndarray:
        areas = np.asarray(self.abs_boundary_area, dtype=float)
        if areas.size == 0 or np.nansum(areas) == 0:
            return np.array([np.nan, np.nan, np.nan], dtype=float)
        centroids = []
        for boundary in self.abs_boundaries:
            centroids.append(_polygon_centroid(boundary))
        centroids = np.asarray(centroids, dtype=float)
        xy = np.nansum(centroids * areas[:, None], axis=0) / np.nansum(areas)
        z = float(np.nansum(self.abs_zPos * areas) / np.nansum(areas))
        return np.array([xy[0], xy[1], z], dtype=float)

    calculate_centroid = CalculateCentroid

    def Plot(self, zInd: int = 1, axisHandle=None):
        import matplotlib.pyplot as plt

        ax = axisHandle if axisHandle is not None else plt.subplots()[1]
        z0 = int(zInd) - 1 if int(zInd) >= 1 else int(zInd)
        handles = []
        if 0 <= z0 < len(self.abs_boundaries):
            b = self.abs_boundaries[z0]
            if len(b):
                (line,) = ax.plot(b[:, 0], b[:, 1])
                handles.append(line)
        ax.set_aspect("equal", adjustable="box")
        return handles

    plot = Plot

    def CalculateMorphology(self):
        eccentricity = np.zeros(self.num_zPos, dtype=float)
        hull_ratio = np.zeros(self.num_zPos, dtype=float)
        num_regions = np.zeros(self.num_zPos, dtype=int)
        for z, boundary in enumerate(self.abs_boundaries):
            if len(boundary) < 3:
                eccentricity[z] = np.nan
                hull_ratio[z] = np.nan
                num_regions[z] = 0
                continue
            cov = np.cov(boundary[:, :2].T)
            eig = np.sort(np.linalg.eigvalsh(cov))[::-1]
            eccentricity[z] = float(math.sqrt(max(0.0, 1.0 - eig[1] / eig[0]))) if eig[0] > 0 else 0.0
            area = _polygon_area(boundary)
            if Polygon is not None:
                poly = _poly(boundary)
                hull_area = float(poly.convex_hull.area) if poly is not None else 0.0
            else:
                hull_area = area
            hull_ratio[z] = float(area / hull_area) if hull_area > 0 else np.nan
            num_regions[z] = 1
        return eccentricity, hull_ratio, num_regions

    calculate_morphology = CalculateMorphology

    def _z_index(self, z_value: float) -> int:
        if self.abs_zPos.size == 0:
            return 0
        candidates = np.where(float(z_value) >= self.abs_zPos)[0]
        if candidates.size == 0:
            return 0
        return int(min(candidates[-1], len(self.abs_boundaries) - 1))

    def DistanceToFeature(self, position: Sequence[float]) -> float:
        pos = np.asarray(position, dtype=float)
        if pos.shape != (3,):
            raise ValueError("position must be a length-3 vector [x, y, z]")
        z = self._z_index(pos[2])
        boundary = self.abs_boundaries[z]
        if len(boundary) == 0:
            return float("inf")
        if Polygon is not None:
            poly = _poly(boundary)
            if poly is not None:
                return float(poly.boundary.distance(Point(float(pos[0]), float(pos[1]))))
        return float(np.min(np.linalg.norm(boundary[:, :2] - pos[:2], axis=1)))

    distance_to_feature = DistanceToFeature

    def IsInFeature(self, position: Sequence[float]) -> bool:
        pos = np.asarray(position, dtype=float)
        if pos.shape != (3,):
            raise ValueError("position must be a length-3 vector [x, y, z]")
        z = self._z_index(pos[2])
        boundary = self.abs_boundaries[z]
        if len(boundary) < 3:
            return False
        if Polygon is not None:
            poly = _poly(boundary)
            return bool(poly is not None and (poly.contains(Point(pos[0], pos[1])) or poly.touches(Point(pos[0], pos[1]))))
        from matplotlib.path import Path as MplPath
        return bool(MplPath(boundary).contains_point(pos[:2]))

    is_in_feature = IsInFeature

    def AssignFeatureID(self, featureID: int) -> None:
        if int(featureID) <= 0 or int(featureID) != featureID:
            raise ValueError("featureID must be a positive integer")
        self.feature_id = int(featureID)

    assign_feature_id = AssignFeatureID

    def DoesFeatureOverlap(self, fFeature: "FoundFeature") -> bool:
        if not isinstance(fFeature, FoundFeature):
            raise TypeError("fFeature must be a FoundFeature")
        for b1, b2 in zip(self.abs_boundaries, fFeature.abs_boundaries):
            if len(b1) < 3 or len(b2) < 3:
                continue
            if Polygon is not None:
                p1 = _poly(b1)
                p2 = _poly(b2)
                if p1 is not None and p2 is not None and p1.intersects(p2):
                    return True
            else:
                if np.any([self.IsInFeature([x, y, self.abs_zPos[0]]) for x, y in b2]):
                    return True
        return False

    does_feature_overlap = DoesFeatureOverlap

    def ReturnPrimaryFovID(self):
        return self.fovID[0] if self.fovID else None

    return_primary_fov_id = ReturnPrimaryFovID

    def CalculateProperties(self) -> None:
        n = len(self.abs_boundaries)
        self.num_zPos = n
        self.boundary_area = np.zeros(n, dtype=float)
        self.abs_boundary_area = np.zeros(n, dtype=float)
        z_thickness = float(np.mean(np.diff(self.abs_zPos))) if len(self.abs_zPos) > 1 else 1.0
        for z in range(n):
            local = self.boundaries[0][z] if self.boundaries and self.boundaries[0] and z < len(self.boundaries[0]) else np.empty((0, 2))
            self.boundary_area[z] = _polygon_area(local)
            self.abs_boundary_area[z] = _polygon_area(self.abs_boundaries[z])
        self.volume = float(np.nansum(self.boundary_area) * z_thickness)
        self.abs_volume = float(np.nansum(self.abs_boundary_area) * z_thickness)

    calculate_properties = CalculateProperties

    def copy(self) -> "FoundFeature":
        cp = _copy.deepcopy(self)
        cp.uID = str(uuid.uuid4())
        return cp

    def CalculateJoinPenalty(self, featureToJoin: Optional["FoundFeature"] = None) -> float:
        if featureToJoin is not None and not isinstance(featureToJoin, FoundFeature):
            raise TypeError("featureToJoin must be a FoundFeature")
        penalties = []
        if featureToJoin is None:
            for b in self.abs_boundaries:
                if len(b):
                    penalties.append(float(np.linalg.norm(b[-1] - b[0])))
        else:
            for b1, b2 in zip(self.abs_boundaries, featureToJoin.abs_boundaries):
                if len(b1) == 0 or len(b2) == 0:
                    continue
                des = np.linalg.norm(b1[-1] - b2[0]) + np.linalg.norm(b1[0] - b2[-1])
                dee = np.linalg.norm(b1[-1] - b2[-1]) + np.linalg.norm(b1[0] - b2[0])
                penalties.append(float(min(des, dee)))
        return float(np.mean(penalties)) if penalties else float("inf")

    calculate_join_penalty = CalculateJoinPenalty

    def JoinFeature(self, featureToJoin: Optional["FoundFeature"] = None) -> "FoundFeature":
        joined = self.copy()
        step = max(float(self.pixel_size) / 1000.0, 1e-6)
        if featureToJoin is None:
            joined.joinedUIDs = [self.uID]
            new_abs = []
            for b in joined.abs_boundaries:
                if len(b) == 0:
                    new_abs.append(b)
                    continue
                gap = _interp_line(b[-1], b[0], step)
                new_abs.append(np.vstack([b, gap]) if len(gap) else b)
            joined.abs_boundaries = new_abs
        else:
            joined.joinedUIDs = [self.uID, featureToJoin.uID]
            joined.fovID = list(self.fovID) + list(featureToJoin.fovID)
            joined.boundaries = [list(self.boundaries[0] if self.boundaries else [])]
            new_abs = []
            for b1, b2 in zip(self.abs_boundaries, featureToJoin.abs_boundaries):
                if len(b1) == 0:
                    new_abs.append(b2)
                    continue
                if len(b2) == 0:
                    new_abs.append(b1)
                    continue
                des = np.linalg.norm(b1[-1] - b2[0]) + np.linalg.norm(b1[0] - b2[-1])
                dee = np.linalg.norm(b1[-1] - b2[-1]) + np.linalg.norm(b1[0] - b2[0])
                if dee < des:
                    b2 = b2[::-1]
                gap12 = _interp_line(b1[-1], b2[0], step)
                gap21 = _interp_line(b2[-1], b1[0], step)
                new_abs.append(np.vstack([b1, gap12, b2, gap21]))
            joined.abs_boundaries = new_abs
        joined.num_joined_features = len(joined.fovID)
        joined.is_broken = False
        joined.CalculateProperties()
        return joined

    join_feature = JoinFeature

    def InFov(self, fovIDs: Iterable[Any]) -> bool:
        fov_set = set(fovIDs if isinstance(fovIDs, (list, tuple, set, np.ndarray)) else [fovIDs])
        return any(f in fov_set for f in self.fovID)

    in_fov = InFov

    def DilateBoundary(self, zIndex: int, dilationSize: float) -> np.ndarray:
        z = int(zIndex) - 1 if int(zIndex) >= 1 else int(zIndex)
        boundary = self.abs_boundaries[z]
        if len(boundary) < 3:
            return boundary.copy()
        if Polygon is not None:
            poly = _poly(boundary)
            if poly is not None:
                dilated = poly.buffer(float(dilationSize))
                if isinstance(dilated, MultiPolygon):
                    dilated = max(dilated.geoms, key=lambda g: g.area)
                return np.asarray(dilated.exterior.coords[:-1], dtype=float)
        centroid = _polygon_centroid(boundary)
        vectors = boundary - centroid
        norms = np.linalg.norm(vectors, axis=1)
        norms[norms == 0] = 1.0
        return boundary + float(dilationSize) * vectors / norms[:, None]

    dilate_boundary = DilateBoundary

    def GeneratePixelMask(self, fovID: Any, zIndices: Sequence[int]) -> np.ndarray:
        if fovID not in self.fovID:
            raise ValueError("requested fovID is not in this feature")
        z_indices = [int(z) - 1 if int(z) >= 1 else int(z) for z in zIndices]
        h, w = self.image_size
        mask = np.zeros((h, w, len(z_indices)), dtype=bool)
        local_idx = self.fovID.index(fovID)
        for out_i, z in enumerate(z_indices):
            b = self.boundaries[local_idx][z] if local_idx < len(self.boundaries) else np.empty((0, 2))
            if len(b) < 3:
                continue
            rr, cc = draw.polygon(b[:, 1], b[:, 0], shape=(h, w))
            mask[rr, cc, out_i] = True
            mask[:, :, out_i] = ndi.binary_fill_holes(mask[:, :, out_i])
        return mask

    generate_pixel_mask = GeneratePixelMask

    def Feature2Table(self) -> pd.DataFrame:
        row = {
            "feature_uID": self.uID,
            "feature_ID": self.feature_id,
            "fovID": self.ReturnPrimaryFovID(),
            "is_broken": self.is_broken,
            "num_joined_features": self.num_joined_features,
            "abs_volume": self.abs_volume,
        }
        for i, boundary in enumerate(self.abs_boundaries, start=1):
            row[f"abs_x_boundary_{i}"] = ";".join(str(float(v)) for v in boundary[:, 0]) if len(boundary) else ""
            row[f"abs_y_boundary_{i}"] = ";".join(str(float(v)) for v in boundary[:, 1]) if len(boundary) else ""
        return pd.DataFrame([row])

    feature2table = Feature2Table
