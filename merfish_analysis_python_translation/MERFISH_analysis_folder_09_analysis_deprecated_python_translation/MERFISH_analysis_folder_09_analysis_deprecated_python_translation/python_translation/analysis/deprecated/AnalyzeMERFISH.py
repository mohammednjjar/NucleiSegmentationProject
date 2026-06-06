
from __future__ import annotations

from pathlib import Path

from .AlignFiducials import align_fiducials
from .CreateWords import create_words
from .DecodeWords import decode_words
from .TransformImageData import transform_image_data
from .helpers import build_image_data_structures, get_field, parse_parameters, read_fasta, read_master_molecule_list, simple_codebook_to_map, transfer_info_file_fields, uuid_string


def analyze_merfish(dataPath, *, parameters=None, imageTag="STORM", imageMListType="alist", fiducialMListType="list", fileExt="bin", fieldNames=None, fieldConv=None, appendExtraFields=True, cellsToAnalyze=None, numHybs=16, bitOrder=None, maxD=8, fiducialFrame=1, fiducialWarp2Hyb1=False, maxDtoCentroid=1, codebookPath="", codebook=None, exactMap=None, correctableMap=None, errCorrFunc=None, keyType="binStr", savePath="", reportsToGenerate=None, verbose=False, printedUpdates=True, image_data_builder=None, molecule_list_reader=None, info_field_transfer=None, composite_image_reporter=None, on_bit_histogram_reporter=None, codebook_mapper=None):
    """Python translation of `AnalyzeMERFISH.m`.

    Runs the legacy MERFISH analysis sequence: discover image/fiducial molecule-list files, load molecule lists, align fiducials, transform image data, construct words, decode with a codebook, and return words plus image/fiducial metadata.
    """
    data_path = Path(dataPath)
    if not data_path.exists() or not data_path.is_dir():
        raise ValueError("A valid data path is required")
    defaults = {
        "imageTag": imageTag,
        "imageMListType": imageMListType,
        "fiducialMListType": fiducialMListType,
        "fileExt": fileExt,
        "fieldNames": fieldNames or ["movieType", "hybNum", "cellNum", "isFiducial", "binType"],
        "fieldConv": fieldConv,
        "appendExtraFields": appendExtraFields,
        "cellsToAnalyze": cellsToAnalyze or [],
        "numHybs": numHybs,
        "bitOrder": bitOrder or list(range(1, int(numHybs) + 1)),
        "maxD": maxD,
        "fiducialFrame": fiducialFrame,
        "fiducialWarp2Hyb1": fiducialWarp2Hyb1,
        "maxDtoCentroid": maxDtoCentroid,
        "codebookPath": codebookPath,
        "codebook": codebook or [],
        "exactMap": exactMap or {},
        "correctableMap": correctableMap or {},
        "errCorrFunc": errCorrFunc,
        "keyType": keyType,
        "savePath": savePath,
        "reportsToGenerate": reportsToGenerate or [],
        "verbose": verbose,
        "printedUpdates": printedUpdates,
    }
    params = parse_parameters(defaults, parameters)
    if params.get("printedUpdates"):
        print("--------------------------------------------------------------")
        print("Analyzing Multiplexed FISH Data")
        print("--------------------------------------------------------------")
        print("Analysis parameters")
        for name in ["codebookPath", "imageTag", "imageMListType", "fiducialMListType", "fileExt", "numHybs", "bitOrder", "maxD", "savePath"]:
            print(" ", name + ":", params.get(name))
    builder = image_data_builder or build_image_data_structures
    reader = molecule_list_reader or read_master_molecule_list
    transfer = info_field_transfer or transfer_info_file_fields
    mapper = codebook_mapper or simple_codebook_to_map
    found_files = builder(str(data_path), parameters=params)
    num_cells = max([int(get_field(f, "cellNum", 0)) for f in found_files], default=0)
    if params.get("printedUpdates"):
        print("--------------------------------------------------------------")
        print("Finding cells in", data_path)
        print(" Found", num_cells, "cells")
    if num_cells == 0:
        raise ValueError("No valid cells found")
    if not params["codebook"] and params.get("codebookPath"):
        params["codebook"] = read_fasta(params["codebookPath"])
    if params["codebook"]:
        params["exactMap"] = mapper(params["codebook"], key_type=params["keyType"])
        if params.get("correctableMap") is None:
            params["correctableMap"] = {}
    if not params.get("exactMap") and params.get("printedUpdates"):
        print("No codebook provided. Found words will not be decoded.")
    requested = params.get("cellsToAnalyze") or list(range(1, num_cells + 1))
    cell_ids = [int(x) for x in requested if 1 <= int(x) <= num_cells]
    words = []
    total_image_data = []
    total_fiducial_data = []
    for cell_id in cell_ids:
        if params.get("printedUpdates"):
            print("--------------------------------------------------------------")
            print("Analyzing data for cell", cell_id, "of", num_cells)
        image_data = [f for f in found_files if get_field(f, "movieType") == params["imageTag"] and get_field(f, "binType") == params["imageMListType"] and int(get_field(f, "cellNum", -1)) == cell_id and not bool(get_field(f, "isFiducial", False))]
        fiducial_data = [f for f in found_files if get_field(f, "movieType") == params["imageTag"] and get_field(f, "binType") == params["fiducialMListType"] and int(get_field(f, "cellNum", -1)) == cell_id and bool(get_field(f, "isFiducial", False))]
        image_data = sorted(image_data, key=lambda x: int(get_field(x, "hybNum", 0)))
        fiducial_data = sorted(fiducial_data, key=lambda x: int(get_field(x, "hybNum", 0)))
        if len(image_data) != len(fiducial_data) or len(image_data) < int(params["numHybs"]) or len(fiducial_data) < int(params["numHybs"]):
            continue
        image_data = image_data[: int(params["numHybs"])]
        fiducial_data = fiducial_data[: int(params["numHybs"])]
        for rec in image_data:
            rec["uID"] = uuid_string()
        for rec in fiducial_data:
            rec["uID"] = uuid_string()
        image_data = transfer(image_data, parameters=params)
        for j in range(int(params["numHybs"])):
            image_data[j]["mList"] = reader(get_field(image_data[j], "filePath"), compact=True, transpose=True, verbose=False)
            fiducial_data[j]["mList"] = reader(get_field(fiducial_data[j], "filePath"), compact=True, transpose=True, verbose=False)
        fiducial_data, params = align_fiducials(fiducial_data, parameters=params)
        image_data, params = transform_image_data(image_data, fiducial_data, parameters=params)
        words_by_cell, params = create_words(image_data, parameters=params)
        if params.get("codebook"):
            words_by_cell, params = decode_words(words_by_cell, params.get("exactMap"), params.get("correctableMap"), parameters=params)
            if params.get("printedUpdates"):
                print(" Found", sum(bool(w.get("isExactMatch")) for w in words_by_cell), "exact matches")
                print(" Found", sum(bool(w.get("isCorrectedMatch")) for w in words_by_cell), "corrected matches")
        if composite_image_reporter is not None:
            composite_image_reporter(words_by_cell, image_data, parameters=params)
        if on_bit_histogram_reporter is not None:
            on_bit_histogram_reporter(words_by_cell, parameters=params)
        words.extend(words_by_cell)
        total_image_data.extend(image_data)
        total_fiducial_data.extend(fiducial_data)
    if params.get("printedUpdates"):
        print("--------------------------------------------------------------")
        print("Completed Multiplexed FISH Analysis")
    return words, total_image_data, total_fiducial_data, params


# MATLAB-style alias
AnalyzeMERFISH = analyze_merfish
