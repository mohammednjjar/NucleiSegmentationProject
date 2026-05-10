"""
Geometric transform utilities for the Python translation of
ZhuangLab/MERFISH_analysis analysis/functions.

The original MATLAB functions use affine2d, fitgeotrans, and
transformPointsForward. This module implements equivalent Python objects for
identity, nonreflective similarity, similarity, affine, projective, and 2-D
polynomial transforms using NumPy least-squares fitting.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple
import math
import numpy as np


ArrayLike = Iterable[Iterable[float]]


def _as_points(points: ArrayLike) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("points must be an N x 2 array or wider; first two columns are used")
    return arr[:, :2]


class GeometricTransform:
    """Base geometric transform interface."""

    kind: str = "identity"

    def transform_points(self, points: ArrayLike) -> np.ndarray:
        return _as_points(points).copy()

    def __call__(self, points: ArrayLike) -> np.ndarray:
        return self.transform_points(points)


@dataclass
class AffineTransform(GeometricTransform):
    """
    2-D affine transform.

    Matrix convention maps row-vector points as:
        [x, y, 1] @ matrix.T -> [x', y']
    where matrix has shape 2 x 3.
    """

    matrix: np.ndarray | None = None
    kind: str = "affine"

    def __post_init__(self) -> None:
        if self.matrix is None:
            self.matrix = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float)
        else:
            self.matrix = np.asarray(self.matrix, dtype=float)
            if self.matrix.shape != (2, 3):
                raise ValueError("affine matrix must have shape 2 x 3")

    def transform_points(self, points: ArrayLike) -> np.ndarray:
        pts = _as_points(points)
        ones = np.ones((pts.shape[0], 1), dtype=float)
        hom = np.hstack([pts, ones])
        return hom @ self.matrix.T


@dataclass
class ProjectiveTransform(GeometricTransform):
    """2-D projective/homography transform for row-vector points."""

    matrix: np.ndarray | None = None
    kind: str = "projective"

    def __post_init__(self) -> None:
        if self.matrix is None:
            self.matrix = np.eye(3, dtype=float)
        else:
            self.matrix = np.asarray(self.matrix, dtype=float)
            if self.matrix.shape != (3, 3):
                raise ValueError("projective matrix must have shape 3 x 3")

    def transform_points(self, points: ArrayLike) -> np.ndarray:
        pts = _as_points(points)
        hom = np.hstack([pts, np.ones((pts.shape[0], 1), dtype=float)])
        mapped = hom @ self.matrix.T
        denom = mapped[:, 2:3]
        with np.errstate(divide="ignore", invalid="ignore"):
            out = mapped[:, :2] / denom
        return out


@dataclass
class PolynomialTransform(GeometricTransform):
    """2-D polynomial transform fitted separately for x' and y'."""

    coefficients_x: np.ndarray
    coefficients_y: np.ndarray
    order: int = 3
    kind: str = "polynomial"

    def _terms(self, points: ArrayLike) -> np.ndarray:
        pts = _as_points(points)
        x = pts[:, 0]
        y = pts[:, 1]
        cols = []
        for total_degree in range(self.order + 1):
            for x_degree in range(total_degree + 1):
                y_degree = total_degree - x_degree
                cols.append((x ** x_degree) * (y ** y_degree))
        return np.vstack(cols).T

    def transform_points(self, points: ArrayLike) -> np.ndarray:
        terms = self._terms(points)
        return np.column_stack([terms @ self.coefficients_x, terms @ self.coefficients_y])


def fit_affine(moving_points: ArrayLike, fixed_points: ArrayLike) -> AffineTransform:
    moving = _as_points(moving_points)
    fixed = _as_points(fixed_points)
    if moving.shape[0] != fixed.shape[0]:
        raise ValueError("moving_points and fixed_points must contain the same number of points")
    if moving.shape[0] < 3:
        raise ValueError("at least 3 control point pairs are required for an affine transform")
    design = np.hstack([moving, np.ones((moving.shape[0], 1), dtype=float)])
    matrix_t, *_ = np.linalg.lstsq(design, fixed, rcond=None)
    return AffineTransform(matrix_t.T)


def fit_nonreflective_similarity(moving_points: ArrayLike, fixed_points: ArrayLike) -> AffineTransform:
    """Fit scale + rotation + translation without reflection."""
    moving = _as_points(moving_points)
    fixed = _as_points(fixed_points)
    if moving.shape[0] != fixed.shape[0]:
        raise ValueError("moving_points and fixed_points must contain the same number of points")
    if moving.shape[0] < 2:
        raise ValueError("at least 2 control point pairs are required for a similarity transform")

    x = moving[:, 0]
    y = moving[:, 1]
    u = fixed[:, 0]
    v = fixed[:, 1]

    zeros = np.zeros_like(x)
    ones = np.ones_like(x)
    design = np.vstack([
        np.column_stack([x, -y, ones, zeros]),
        np.column_stack([y,  x, zeros, ones]),
    ])
    target = np.concatenate([u, v])
    a, b, tx, ty = np.linalg.lstsq(design, target, rcond=None)[0]
    matrix = np.array([[a, -b, tx], [b, a, ty]], dtype=float)
    return AffineTransform(matrix, kind="nonreflectivesimilarity")


def fit_similarity(moving_points: ArrayLike, fixed_points: ArrayLike, allow_reflection: bool = True) -> AffineTransform:
    """
    Fit a similarity transform using Umeyama alignment.

    MATLAB's 'similarity' can represent reflection; this implementation tries
    both non-reflective and reflective solutions when allow_reflection=True and
    returns the lower-RMSE solution.
    """
    moving = _as_points(moving_points)
    fixed = _as_points(fixed_points)
    if moving.shape[0] != fixed.shape[0]:
        raise ValueError("moving_points and fixed_points must contain the same number of points")
    if moving.shape[0] < 2:
        raise ValueError("at least 2 control point pairs are required for a similarity transform")

    candidates = [fit_nonreflective_similarity(moving, fixed)]
    if allow_reflection:
        reflected = moving.copy()
        reflected[:, 0] *= -1.0
        t_ref = fit_nonreflective_similarity(reflected, fixed)
        reflection_matrix = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        affine_3 = np.vstack([t_ref.matrix, [0.0, 0.0, 1.0]]) @ reflection_matrix
        candidates.append(AffineTransform(affine_3[:2, :], kind="similarity"))

    def rmse(transform: AffineTransform) -> float:
        diff = transform.transform_points(moving) - fixed
        return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))

    best = min(candidates, key=rmse)
    best.kind = "similarity" if allow_reflection else "nonreflectivesimilarity"
    return best


def fit_projective(moving_points: ArrayLike, fixed_points: ArrayLike) -> ProjectiveTransform:
    moving = _as_points(moving_points)
    fixed = _as_points(fixed_points)
    if moving.shape[0] != fixed.shape[0]:
        raise ValueError("moving_points and fixed_points must contain the same number of points")
    if moving.shape[0] < 4:
        raise ValueError("at least 4 control point pairs are required for a projective transform")

    rows = []
    for (x, y), (u, v) in zip(moving, fixed):
        rows.append([-x, -y, -1.0, 0.0, 0.0, 0.0, u * x, u * y, u])
        rows.append([0.0, 0.0, 0.0, -x, -y, -1.0, v * x, v * y, v])
    a = np.asarray(rows, dtype=float)
    _, _, vh = np.linalg.svd(a)
    h = vh[-1].reshape(3, 3)
    if abs(h[2, 2]) > np.finfo(float).eps:
        h = h / h[2, 2]
    return ProjectiveTransform(h)


def _polynomial_terms(points: np.ndarray, order: int) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    cols = []
    for total_degree in range(order + 1):
        for x_degree in range(total_degree + 1):
            y_degree = total_degree - x_degree
            cols.append((x ** x_degree) * (y ** y_degree))
    return np.vstack(cols).T


def fit_polynomial(moving_points: ArrayLike, fixed_points: ArrayLike, order: int = 3) -> PolynomialTransform:
    moving = _as_points(moving_points)
    fixed = _as_points(fixed_points)
    if moving.shape[0] != fixed.shape[0]:
        raise ValueError("moving_points and fixed_points must contain the same number of points")
    n_terms = (order + 1) * (order + 2) // 2
    if moving.shape[0] < n_terms:
        raise ValueError(f"at least {n_terms} control point pairs are required for polynomial order {order}")
    terms = _polynomial_terms(moving, order)
    cx, *_ = np.linalg.lstsq(terms, fixed[:, 0], rcond=None)
    cy, *_ = np.linalg.lstsq(terms, fixed[:, 1], rcond=None)
    return PolynomialTransform(cx, cy, order=order)


def fit_geometric_transform(
    moving_points: ArrayLike,
    fixed_points: ArrayLike,
    transformation_type: str = "nonreflectivesimilarity",
    polynomial_order: int = 3,
) -> GeometricTransform:
    """Python equivalent of MATLAB fitgeotrans for the types used here."""
    kind = transformation_type.lower()
    if kind == "nonreflectivesimilarity":
        return fit_nonreflective_similarity(moving_points, fixed_points)
    if kind == "similarity":
        return fit_similarity(moving_points, fixed_points, allow_reflection=True)
    if kind == "affine":
        return fit_affine(moving_points, fixed_points)
    if kind == "projective":
        return fit_projective(moving_points, fixed_points)
    if kind == "polynomial":
        return fit_polynomial(moving_points, fixed_points, order=int(polynomial_order))
    raise ValueError(f"unsupported transformation_type: {transformation_type!r}")


def transform_points_forward(transform: GeometricTransform | None, points: ArrayLike) -> np.ndarray:
    if transform is None:
        return _as_points(points).copy()
    return transform.transform_points(points)


def identity_transform() -> AffineTransform:
    return AffineTransform()
