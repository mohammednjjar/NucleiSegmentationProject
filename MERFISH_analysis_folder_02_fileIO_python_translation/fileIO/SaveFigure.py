from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any
import warnings

from . import _utils
from ._utils import make_parameters


def SaveFigure(figHandle, **kwargs: Any):
    """Save a Matplotlib figure in one or more formats.

    MATLAB .fig output is translated as a pickled Matplotlib figure stored
    with .fig extension. The MATLAB export_fig dependency is replaced by
    Matplotlib's native savefig.
    """
    defaults = {
        "addGit": False,
        "gitState": None,
        "name": None,
        "overwrite": False,
        "verbose": True,
        "formats": ["fig"],
        "appendFormats": [],
        "useExportFig": True,
        "subFolder": None,
        "makeDir": True,
        "savePath": _utils.FIGURE_SAVE_PATH,
        "transparent": False,
        "closeFig": False,
        "saveData": True,
    }
    parameters = make_parameters(defaults, kwargs)
    if figHandle is None or not hasattr(figHandle, "savefig"):
        raise ValueError("Provided figure handle is not valid")

    save_path = Path(parameters.savePath or _utils.FIGURE_SAVE_PATH or ".")
    if parameters.subFolder:
        save_path = save_path / parameters.subFolder
    if parameters.makeDir:
        save_path.mkdir(parents=True, exist_ok=True)
    elif not save_path.is_dir():
        raise ValueError("Provided subfolder path does not exist or cannot be created")

    name = parameters.name or getattr(figHandle, "_suptitle", None)
    if name is not None and not isinstance(name, str):
        name = getattr(name, "get_text", lambda: "unnamed")()
    name = name or "unnamed"

    formats = list(parameters.formats or [])
    for fmt in parameters.appendFormats or []:
        if fmt not in formats:
            formats.append(fmt)

    if not parameters.overwrite:
        count = 0
        new_name = name
        while any((save_path / f"{new_name}.{fmt}").exists() for fmt in formats):
            count += 1
            new_name = f"{name}({count})"
        if count > 0 and parameters.verbose:
            print(f"Found existing files. Appending {count} to all file names.")
        name = new_name

    saved_file_paths = []
    if parameters.saveData:
        for fmt in formats:
            out = save_path / f"{name}.{fmt}"
            if fmt == "fig":
                with out.open("wb") as fh:
                    pickle.dump(figHandle, fh, protocol=pickle.HIGHEST_PROTOCOL)
            elif fmt == "ai":
                figHandle.savefig(out, format="eps", transparent=parameters.transparent)
            elif fmt in {"eps", "png", "pdf", "svg", "jpg", "jpeg", "tif", "tiff"}:
                figHandle.savefig(out, format=fmt, transparent=parameters.transparent)
            else:
                warnings.warn(f"Unrecognized/unsupported format {fmt}", RuntimeWarning)
                continue
            saved_file_paths.append(str(out))
            if parameters.verbose:
                print(f"Saved: {out}")

    if parameters.closeFig:
        try:
            import matplotlib.pyplot as plt
            plt.close(figHandle)
        except Exception as exc:
            _ = exc
    return saved_file_paths, parameters
