"""Path resolution shared by the pipeline, the revision analyses and the figure scripts.

The working copy these analyses were developed in keeps pipeline outputs under
``results/tables/`` with a curated copy under ``paper_handoff/tables/``; the public release
flattens both into ``tables/``. Hard-coding either layout makes the released scripts fail
from a fresh clone, so resolve at run time instead.

Usage::

    from src.paths import ROOT, find_table, table_out_dirs

    df = pd.read_csv(find_table("oof_predictions.tsv"), sep="\\t")
    for d in table_out_dirs():
        df.to_csv(os.path.join(d, "my_output.tsv"), sep="\\t", index=False)
"""
from __future__ import annotations

import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# Searched in order; the first that exists wins.
TABLE_DIRS = ("results/tables", "paper_handoff/tables", "tables")


def find_table(name: str) -> str:
    """Absolute path to a released result table, whichever layout this checkout uses."""
    tried = []
    for sub in TABLE_DIRS:
        p = os.path.join(ROOT, sub, name)
        tried.append(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        "%s not found in this checkout; tried:\n  %s" % (name, "\n  ".join(tried))
    )


def table_out_dirs() -> list:
    """Directories to write pipeline tables into.

    Writes to every table directory that already exists, so the working copy keeps both
    its canonical and curated copies in step. In a fresh clone only ``tables/`` exists,
    and it is created if missing.
    """
    dirs = [os.path.join(ROOT, s) for s in TABLE_DIRS if os.path.isdir(os.path.join(ROOT, s))]
    if not dirs:
        dirs = [os.path.join(ROOT, "tables")]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    return dirs


def revision_out_dir() -> str:
    """Directory the revision analyses write their tables into (alongside the scripts)."""
    d = os.path.join(ROOT, "scripts", "revision")
    os.makedirs(d, exist_ok=True)
    return d
