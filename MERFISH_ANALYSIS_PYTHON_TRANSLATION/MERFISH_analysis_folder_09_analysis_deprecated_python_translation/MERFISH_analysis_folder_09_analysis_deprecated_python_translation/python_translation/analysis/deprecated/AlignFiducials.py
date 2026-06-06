
from __future__ import annotations

import numpy as np

from .Warp2BestPair import warp2_best_pair
from .helpers import AffineTransform, compose_transforms, filter_mlist, get_array_field, get_field, mlist_xy, parse_parameters, set_field


def align_fiducials(fiducialData, *, parameters=None, maxD=8, fiducialFrame=1, fiducialWarp2Hyb1=False, verbose=True, printedUpdates=True, reportsToGenerate=None, useSubFolderForCellReport=True, overwrite=True, figFormats=("png", "fig"), saveAndClose=False):
    """Python translation of `AlignFiducials.m`.

    Builds affine/rigid transforms from fiducial molecule lists. Each input record must contain an `mList` field with x/y coordinates (`xc`, `yc`). The returned records include `tform`, `warpErrors`, `hasFiducialError`, and `fiducialErrorMessage`.
    """
    data = list(fiducialData)
    if len(data) < 2:
        raise ValueError("Invalid fiducial data: at least two fiducial records are required")
    defaults = {
        "maxD": maxD,
        "fiducialFrame": fiducialFrame,
        "fiducialWarp2Hyb1": fiducialWarp2Hyb1,
        "verbose": verbose,
        "printedUpdates": printedUpdates,
        "reportsToGenerate": reportsToGenerate or [],
        "useSubFolderForCellReport": useSubFolderForCellReport,
        "overwrite": overwrite,
        "figFormats": figFormats,
        "saveAndClose": saveAndClose,
        "numHybs": len(data),
    }
    params = parse_parameters(defaults, parameters)
    if params.get("printedUpdates") and params.get("verbose"):
        print("--------------------------------------------------------------")
        print("Analyzing fiducials")
    fed_pos = []
    for rec in data:
        mlist = filter_mlist(get_field(rec, "mList"), frame=[1, 2])
        fed_pos.append(mlist_xy(mlist, "xc", "yc"))
    nearest_transforms: list[AffineTransform | None] = [None] * len(data)
    for i, rec in enumerate(data):
        try:
            set_field(rec, "warpErrors", np.full(5, np.nan))
            if bool(params.get("fiducialWarp2Hyb1")):
                transform, errors = warp2_best_pair(fed_pos[0], fed_pos[i], parameters=params, showPlots=False)
            else:
                reference_i = max(i - 1, 0)
                nearest_transform, errors = warp2_best_pair(fed_pos[reference_i], fed_pos[i], parameters=params, showPlots=False)
                nearest_transforms[i] = nearest_transform
                transform = compose_transforms([t for t in nearest_transforms[: i + 1] if t is not None])
            set_field(rec, "tform", transform)
            set_field(rec, "hasFiducialError", False)
            set_field(rec, "fiducialErrorMessage", "")
            clean = np.full(5, np.nan)
            errors = np.asarray(errors, dtype=float).ravel()
            clean[: min(len(errors), 5)] = errors[: min(len(errors), 5)]
            set_field(rec, "warpErrors", clean)
        except Exception as exc:
            if i == 0:
                transform = AffineTransform.identity()
            else:
                transform = get_field(data[i - 1], "tform", AffineTransform.identity())
            set_field(rec, "tform", transform)
            set_field(rec, "warpErrors", np.full(5, np.nan))
            set_field(rec, "hasFiducialError", True)
            set_field(rec, "fiducialErrorMessage", str(exc))
        if params.get("verbose"):
            print("Alignment error =", get_field(rec, "warpErrors"))
    return data, params


# MATLAB-style alias
AlignFiducials = align_fiducials
