
from __future__ import annotations

from .helpers import get_array_field, get_field, parse_parameters, set_field


def transform_image_data(imageData, fiducialData, *, parameters=None, verbose=True, printedUpdates=True):
    """Python translation of `TransformImageData.m`.

    Applies each fiducial inverse transform to the corresponding image molecule-list `x/y` coordinates and writes aligned `xc/yc` coordinates, then transfers transform/error/UID metadata.
    """
    images = list(imageData)
    fiducials = list(fiducialData)
    if len(images) != len(fiducials):
        raise ValueError("imageData and fiducialData must have the same length")
    params = parse_parameters({"verbose": verbose, "printedUpdates": printedUpdates}, parameters)
    if params.get("printedUpdates") and params.get("verbose"):
        print("--------------------------------------------------------------")
        print("Shifting images")
    fields_from = ["tform", "warpErrors", "hasFiducialError", "fiducialErrorMessage", "uID"]
    fields_to = ["tform", "warpErrors", "hasFiducialError", "fiducialErrorMessage", "fidUID"]
    for img, fid in zip(images, fiducials):
        mlist = get_field(img, "mList")
        x = get_array_field(mlist, "x")
        y = get_array_field(mlist, "y")
        points = get_field(fid, "tform").inverse_points(__import__("numpy").column_stack([x, y]))
        set_field(mlist, "xc", points[:, 0])
        set_field(mlist, "yc", points[:, 1])
        for src, dst in zip(fields_from, fields_to):
            set_field(img, dst, get_field(fid, src))
    return images, params


# MATLAB-style alias
TransformImageData = transform_image_data
