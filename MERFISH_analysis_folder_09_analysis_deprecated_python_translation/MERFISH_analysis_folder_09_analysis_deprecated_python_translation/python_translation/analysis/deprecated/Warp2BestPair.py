
from __future__ import annotations

import numpy as np

from .helpers import AffineTransform, estimate_rigid_transform, mutual_nearest_pairs, nearest_neighbors, parse_parameters


def match_fiducials(hybe1, hybe2, *, parameters=None, maxD=None):
    """Match fiducial bead positions between two hybridization images.

    MATLAB used `MatchFeducials`; this Python translation performs mutual nearest-neighbor matching within `maxD`.
    If too few mutual matches are found, it falls back to nearest-neighbor matches sorted by distance.
    """
    params = parse_parameters({"maxD": 2.0}, parameters)
    if maxD is not None:
        params["maxD"] = maxD
    reference = np.asarray(hybe1, dtype=float).reshape(-1, 2)
    moving = np.asarray(hybe2, dtype=float).reshape(-1, 2)
    ref_ids, mov_ids, dist = mutual_nearest_pairs(reference, moving, float(params["maxD"]))
    if len(ref_ids) >= 2:
        return ref_ids, mov_ids, params
    idx_ref, nn_dist = nearest_neighbors(moving, reference)
    order = np.argsort(nn_dist)
    ref_ids = []
    mov_ids = []
    used = set()
    for moving_i in order:
        ref_i = int(idx_ref[moving_i])
        if ref_i < 0 or ref_i in used:
            continue
        used.add(ref_i)
        ref_ids.append(ref_i)
        mov_ids.append(int(moving_i))
        if len(ref_ids) >= max(2, min(len(reference), len(moving))):
            break
    return np.asarray(ref_ids, dtype=int), np.asarray(mov_ids, dtype=int), params


def warp2_best_pair(hybe1, hybe2, *, parameters=None, maxD=2.0, useCorrAlign=True, fighandle=None, imageSize=(256, 256), showPlots=True, verbose=True):
    """Python translation of `Warp2BestPair.m`.

    Computes a rotation+translation transform by matching fiducial points in two images,
    selecting the pair of matches with the most mutually consistent shift vectors, and
    estimating a rigid transform. The returned transform maps reference (`hybe1`) points
    to moving (`hybe2`) points; `transform.inverse_points(hybe2)` aligns moving points
    back into the reference frame, matching MATLAB `tforminv` behavior.
    """
    params = parse_parameters({
        "maxD": maxD,
        "useCorrAlign": useCorrAlign,
        "fighandle": fighandle,
        "imageSize": imageSize,
        "showPlots": showPlots,
        "verbose": verbose,
    }, parameters)
    reference = np.asarray(hybe1, dtype=float).reshape(-1, 2)
    moving = np.asarray(hybe2, dtype=float).reshape(-1, 2)
    if reference.shape[1] != 2 or moving.shape[1] != 2:
        raise ValueError("hybe1 and hybe2 must be n x 2 point arrays")
    matched1, matched2, params = match_fiducials(reference, moving, parameters=params)
    if len(matched1) < 2:
        raise ValueError("found fewer than 2 fiducials, cannot compute warp")
    shifts = reference[matched1] - moving[matched2]
    if len(shifts) > 2:
        diff = shifts[:, None, :] - shifts[None, :, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2))
        np.fill_diagonal(dist, np.inf)
        row, col = np.unravel_index(np.argmin(dist), dist.shape)
        pair_ref = np.asarray([matched1[row], matched1[col]], dtype=int)
        pair_mov = np.asarray([matched2[row], matched2[col]], dtype=int)
    else:
        pair_ref = matched1[:2]
        pair_mov = matched2[:2]
    transform = estimate_rigid_transform(reference[pair_ref], moving[pair_mov])
    warped_moving = transform.inverse_points(moving)
    _, distances = nearest_neighbors(warped_moving, reference)
    distances = np.asarray(distances, dtype=float)
    distances[distances > float(params["maxD"])] = np.nan
    warp_errors = np.sort(distances)
    return transform, warp_errors


# MATLAB-style alias
Warp2BestPair = warp2_best_pair
