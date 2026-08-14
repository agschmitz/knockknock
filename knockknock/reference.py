"""Reference data loaders (Matveeva coefficients, RNA NN ΔG, ΔG bins).

All reference tables ship as CSV in ``knockknock/data/`` so they remain
human-readable and immune to pickle version drift.  See ``data/README.md``
for source publications and column meanings.
"""

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path

import pandas as pd


def get_data_dir() -> Path:
    return Path(__file__).parent / "data"


def load_coeffs() -> pd.DataFrame:
    """Position-specific nucleotide efficacy coefficients (Matveeva 2012).

    Index: ``position`` (0-based). Columns: ``A, C, G, T, dG``. The engine only
    consumes the four nucleotide columns; ``dG`` is preserved for completeness.
    """
    return pd.read_csv(get_data_dir() / "coeffs.csv", index_col="position")


def load_deltaG() -> dict[str, float]:
    """RNA nearest-neighbor ΔG (Xia et al. 1998).

    The on-disk table only stores one orientation of each NN stack. We expand
    to all 16 RNA dinucleotides using nearest-neighbor symmetry: the stack
    5'-XY-3' has the same ΔG as its reverse-complement 5'-rc(Y)+rc(X)-3'.
    """
    df = pd.read_csv(get_data_dir() / "deltaG.csv")
    full: dict[str, float] = dict(
        zip(df["dinucleotide"], df["dG"].astype(float), strict=False)
    )

    rc_map = str.maketrans("ACGU", "UGCA")
    for x, y in product("ACGU", repeat=2):
        key = x + y
        if key in full:
            continue
        rc = (y + x).translate(rc_map)
        if rc in full:
            full[key] = full[rc]
    if len(full) != 16:
        missing = sorted({a + b for a in "ACGU" for b in "ACGU"} - set(full.keys()))
        raise RuntimeError(
            f"deltaG table incomplete after symmetry expansion: missing {missing}"
        )
    return full


def load_dG_bins() -> pd.DataFrame:
    """Position-binned ΔG penalty table.

    Index: ``position`` (0-based). Columns: ``bin_1..bin_4`` corresponding to
    the four ΔG-NN stability classes used in :func:`engine._score_bins`.
    """
    return pd.read_csv(get_data_dir() / "dG_bins.csv", index_col="position")


@dataclass
class ReferenceData:
    """Bundle of all reference tables consumed by the scoring engine."""

    coeffs: pd.DataFrame
    deltaG: dict[str, float]
    dG_bins: pd.DataFrame


@lru_cache(maxsize=1)
def load_reference_data() -> ReferenceData:
    """Load all three reference tables once, cached for the process lifetime."""
    return ReferenceData(
        coeffs=load_coeffs(),
        deltaG=load_deltaG(),
        dG_bins=load_dG_bins(),
    )
