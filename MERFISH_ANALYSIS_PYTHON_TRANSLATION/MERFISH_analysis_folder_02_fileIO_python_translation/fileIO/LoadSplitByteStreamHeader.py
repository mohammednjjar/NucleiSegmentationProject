from __future__ import annotations

from .LoadByteStream import LoadByteStream


def LoadSplitByteStreamHeader(filePath: str):
    """Load and validate split-byte-stream metadata."""
    sbs_info, parameters = LoadByteStream(filePath, verbose=False)
    if not isinstance(sbs_info, dict) or sbs_info.get("isByteStream") != "Y":
        raise ValueError("The specified matb file does not represent a split byte stream")
    return sbs_info, parameters
