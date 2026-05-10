from __future__ import annotations

from typing import Any

from ..BuildFileStructure import BuildFileStructure
from .CreateImageDataStructure import CreateImageDataStructure


def BuildImageDataStructures(folderPath: str, **kwargs: Any):
    """Build imageData dictionaries for files matching BuildFileStructure rules."""
    found_files, parameters = BuildFileStructure(folderPath, **kwargs)
    default = CreateImageDataStructure(1)[0]
    image_data = []
    for record in found_files:
        merged = dict(default)
        merged.update(record)
        image_data.append(merged)
    return image_data, parameters
