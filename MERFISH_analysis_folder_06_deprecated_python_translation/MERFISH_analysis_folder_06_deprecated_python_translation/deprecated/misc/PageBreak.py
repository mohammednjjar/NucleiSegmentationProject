"""Python translation of `deprecated/misc/PageBreak.m`."""
from __future__ import annotations


def PageBreak(display: bool = True) -> str:
    break_string = '-------------------------------------------------------------------------'
    if display:
        print(break_string)
    return break_string


def page_break(display: bool = True) -> str:
    return PageBreak(display=display)
