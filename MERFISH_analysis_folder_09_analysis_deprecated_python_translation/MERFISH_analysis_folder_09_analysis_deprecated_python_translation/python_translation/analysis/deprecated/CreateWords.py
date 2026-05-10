
from __future__ import annotations

import numpy as np

from .helpers import create_words_structure, get_array_field, get_field, has_field, mlist_xy, nearest_neighbors, parse_parameters, set_field, stable_unique_rows, uuid_string, bits_to_int


def _connected_components(mask: np.ndarray) -> list[np.ndarray]:
    visited = np.zeros(mask.shape, dtype=bool)
    comps: list[np.ndarray] = []
    height, width = mask.shape
    for r in range(height):
        for c in range(width):
            if not mask[r, c] or visited[r, c]:
                continue
            stack = [(r, c)]
            visited[r, c] = True
            coords = []
            while stack:
                cr, cc = stack.pop()
                coords.append((cr, cc))
                for nr in (cr - 1, cr, cr + 1):
                    for nc in (cc - 1, cc, cc + 1):
                        if nr == cr and nc == cc:
                            continue
                        if 0 <= nr < height and 0 <= nc < width and mask[nr, nc] and not visited[nr, nc]:
                            visited[nr, nc] = True
                            stack.append((nr, nc))
            comps.append(np.asarray(coords, dtype=int))
    return comps


def _copy_word_metadata(word, imageData, mLists, mListFields, params):
    bit_order = np.asarray(params.get("bitOrder", np.arange(1, int(params["numHybs"]) + 1)), dtype=int) - 1
    bit_order = bit_order[(bit_order >= 0) & (bit_order < int(params["numHybs"]))]
    measured = np.asarray(word["measuredCodeword"], dtype=bool).ravel()
    codeword = measured[bit_order]
    word["imageNames"] = [get_field(img, "name", "") for img in imageData]
    word["imagePaths"] = [get_field(img, "filePath", "") for img in imageData]
    word["bitOrder"] = (bit_order + 1).tolist()
    word["numHyb"] = int(params["numHybs"])
    word["cellID"] = get_field(imageData[0], "cellNum", 1) if imageData else 1
    word["imageX"] = get_field(imageData[0], "Stage_X", 0.0) if imageData else 0.0
    word["imageY"] = get_field(imageData[0], "Stage_Y", 0.0) if imageData else 0.0
    word["uID"] = uuid_string()
    word["imageUIDs"] = [get_field(img, "uID", "") for img in imageData]
    word["fiducialUIDs"] = [get_field(img, "fidUID", "") for img in imageData]
    word["hasFiducialError"] = [bool(get_field(img, "hasFiducialError", False)) for img in imageData]
    word["codeword"] = codeword
    word["intCodeword"] = bits_to_int(codeword)
    word["numOnBits"] = int(np.sum(measured))
    word["measuredOnBits"] = (np.flatnonzero(measured) + 1).tolist()
    word["onBits"] = (np.flatnonzero(codeword) + 1).tolist()
    word["paddedCellID"] = [int(word["cellID"])] * int(word["numOnBits"])
    for list_index, molecule_index in zip(np.flatnonzero(measured), list(word.get("mListInds", []))):
        if np.isnan(molecule_index):
            continue
        molecule_i = int(molecule_index)
        mlist = mLists[list_index]
        for field in mListFields:
            values = get_field(mlist, field, [])
            arr = np.asarray(values)
            if arr.ndim >= 1 and molecule_i < len(arr):
                target = word.setdefault(field, [np.nan] * len(mLists))
                if not isinstance(target, list):
                    target = [np.nan] * len(mLists)
                    word[field] = target
                target[list_index] = arr[molecule_i].item() if hasattr(arr[molecule_i], "item") else arr[molecule_i]
    word["geneName"] = ""
    word["isExactMatch"] = False
    word["isCorrectedMatch"] = False
    word["focusLockQuality"] = [get_field(img, "focusLockQuality", np.nan) for img in imageData]


def create_words(imageData, *, parameters=None, wordConstMethod="perLocalization", binSize=0.25, minDotPerBin=1, minLocsPerDot=1, minArea=0, maxArea=10, showPlots=False, verbose=True, clusterFig=None, histFig=None, imageSize=(256, 256), printedUpdates=True, numHybs=16, reportsToGenerate=None, useSubFolderForCellReport=True, overwrite=True, figFormats=("png", "fig"), maxDtoCentroid=1, bitOrder=None):
    """Python translation of `CreateWords.m`.

    Reconstructs MERFISH words from aligned molecule lists using either `commonCentroid`
    connected-component binning or `perLocalization` nearest-neighbor word construction.
    """
    images = list(imageData)
    defaults = {
        "wordConstMethod": wordConstMethod,
        "binSize": binSize,
        "minDotPerBin": minDotPerBin,
        "minLocsPerDot": minLocsPerDot,
        "minArea": minArea,
        "maxArea": maxArea,
        "showPlots": showPlots,
        "verbose": verbose,
        "clusterFig": clusterFig,
        "histFig": histFig,
        "imageSize": imageSize,
        "printedUpdates": printedUpdates,
        "numHybs": numHybs,
        "reportsToGenerate": reportsToGenerate or [],
        "useSubFolderForCellReport": useSubFolderForCellReport,
        "overwrite": overwrite,
        "figFormats": figFormats,
        "maxDtoCentroid": maxDtoCentroid,
        "bitOrder": bitOrder or list(range(1, int(numHybs) + 1)),
    }
    params = parse_parameters(defaults, parameters)
    if params.get("printedUpdates") and params.get("verbose"):
        print("--------------------------------------------------------------")
        print("Creating words with method:", params["wordConstMethod"])
    mLists = [get_field(img, "mList", {}) for img in images]
    mListFields = list(mLists[0].keys()) if mLists and isinstance(mLists[0], dict) else []
    words = []
    if params["wordConstMethod"] == "commonCentroid":
        point_arrays = [mlist_xy(mlist, "xc", "yc") for mlist in mLists]
        point_arrays = [arr for arr in point_arrays if len(arr)]
        all_points = np.vstack(point_arrays) if point_arrays else np.zeros((0, 2))
        if len(all_points):
            bin_size = float(params["binSize"])
            image_size = np.asarray(params["imageSize"], dtype=float).ravel()
            y_edges = np.arange(bin_size, image_size[0] + bin_size, bin_size)
            x_edges = np.arange(bin_size, image_size[1] + bin_size, bin_size)
            hist, _y_edges, _x_edges = np.histogram2d(all_points[:, 1], all_points[:, 0], bins=[y_edges, x_edges])
            components = _connected_components(hist >= float(params["minDotPerBin"]))
            centroids = []
            for comp in components:
                values = hist[comp[:, 0], comp[:, 1]]
                area = len(comp)
                locs = float(np.sum(values))
                if locs > float(params["minLocsPerDot"]) and area > float(params["minArea"]) and area < float(params["maxArea"]):
                    weighted_yx = np.average(comp.astype(float), axis=0, weights=values)
                    centroids.append([weighted_yx[1] * bin_size, weighted_yx[0] * bin_size])
            putative = np.asarray(centroids, dtype=float).reshape(-1, 2)
        else:
            putative = np.zeros((0, 2))
        idx = np.zeros((len(putative), len(images)), dtype=int)
        dist = np.full((len(putative), len(images)), np.inf)
        for i, mlist in enumerate(mLists):
            pts = mlist_xy(mlist, "xc", "yc")
            if len(pts) and len(putative):
                idx[:, i], dist[:, i] = nearest_neighbors(putative, pts)
        words = create_words_structure(len(putative), int(params["numHybs"]))
        for i, centroid in enumerate(putative):
            measured = dist[i, :] <= float(params["maxDtoCentroid"])
            words[i]["measuredCodeword"] = measured
            words[i]["mListInds"] = idx[i, measured].astype(float).tolist()
            words[i]["wordCentroidX"] = float(centroid[0])
            words[i]["wordCentroidY"] = float(centroid[1])
    elif params["wordConstMethod"] == "perLocalization":
        min_phots_per_stain = 1
        spot_positions = [mlist_xy(mlist, "xc", "yc") for mlist in mLists]
        if spot_positions and all(len(x) > 0 for x in spot_positions):
            putative = np.vstack(spot_positions)
            num_hybes = len(images)
            detected = np.zeros((len(putative), num_hybes), dtype=bool)
            x_per = np.full((len(putative), num_hybes), np.nan)
            y_per = np.full((len(putative), num_hybes), np.nan)
            idx_per = np.full((len(putative), num_hybes), np.nan)
            for h, mlist in enumerate(mLists):
                pts = spot_positions[h]
                idx, di = nearest_neighbors(putative, pts)
                brightness = get_array_field(mlist, "a", np.ones(len(pts)))
                if len(brightness) < len(pts):
                    brightness = np.pad(brightness, (0, len(pts) - len(brightness)), constant_values=1)
                b = brightness[idx]
                b[di >= float(params["maxDtoCentroid"])] = 0
                valid = b > min_phots_per_stain
                detected[:, h] = valid
                x_per[:, h] = pts[idx, 0]
                y_per[:, h] = pts[idx, 1]
                x_per[~valid, h] = np.nan
                y_per[~valid, h] = np.nan
                idx_per[:, h] = idx.astype(float)
                idx_per[~valid, h] = np.nan
            mean_x = np.nanmean(x_per, axis=1)
            mean_y = np.nanmean(y_per, axis=1)
            locations, unique_idx = stable_unique_rows(np.column_stack([mean_x, mean_y]))
            detected = detected[unique_idx]
            idx_per = idx_per[unique_idx]
            words = create_words_structure(detected.shape[0], int(params["numHybs"]))
            for i in range(detected.shape[0]):
                words[i]["measuredCodeword"] = detected[i, :]
                words[i]["mListInds"] = idx_per[i, ~np.isnan(idx_per[i, :])].astype(float).tolist()
                words[i]["wordCentroidX"] = float(locations[i, 0])
                words[i]["wordCentroidY"] = float(locations[i, 1])
        else:
            words = create_words_structure(0, int(params["numHybs"]))
    else:
        raise ValueError("Unknown word construction method")
    for i, word in enumerate(words):
        word["wordNumInCell"] = i + 1
        _copy_word_metadata(word, images, mLists, mListFields, params)
    if params.get("printedUpdates"):
        print(" Reconstructed", len(words), "words")
    return words, params


# MATLAB-style alias
CreateWords = create_words
