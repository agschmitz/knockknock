"""Backbone definitions — dataclass and CSV loader.

Junctions are stored as `cutsite.Overhang` objects derived from the backbone's
cut sites, written in conventional notation in ``data/backbones.csv``
(``G/AATTC``, ``CTGCA/G``, ...). See `knockknock.cutsite` for the notation and
the symmetric-cut rule that turns it into an overhang.

The ``5_prime_site``/``3_prime_site`` columns describe the **digested
backbone**, never the insert — which is why the two columns are parsed with
their own `LEFT`/`RIGHT` junction. It only changes the reading of a directly
stated overhang (``5'TTAT``), and only for non-palindromic ones, but getting
it wrong there designs an oligo that cannot ligate.
"""

import csv
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from knockknock.cutsite import LEFT, RIGHT, CutSiteError, Overhang, parse_overhang


@dataclass
class Backbone:
    name: str
    five_prime_enzyme: str
    #: The backbone's left-hand end, i.e. the one the insert's 5' (left)
    #: junction has to fit.
    five_prime_overhang: Overhang
    five_prime_flank: str
    loop: str
    three_prime_flank: str
    three_prime_enzyme: str
    #: The backbone's right-hand end, i.e. the one the insert's 3' (right)
    #: junction has to fit.
    three_prime_overhang: Overhang
    guide_on_5p: bool
    #: Carry the rest of each recognition sequence on the insert, so ligation
    #: rebuilds the restriction site instead of scarring it. A property of how
    #: this vector is normally cloned into, so it travels with the backbone —
    #: `--keep-sites`/`--no-keep-sites` still override it per run. Defaults to
    #: False for a row whose CSV predates the column.
    keep_sites: bool = False
    #: 1-based passenger-strand mismatch position(s), counted from the base
    #: of the stem (see `oligo.compute_mismatched_passenger`). A property of
    #: how this vector is conventionally opened for Drosha/Dicer processing,
    #: so it travels with the backbone the same way `keep_sites` does —
    #: `-m`/`--mismatch` still overrides it per run. Empty means "no
    #: mismatch" (a perfect-complement passenger). Defaults to ``[1]`` for a
    #: row whose CSV predates the column.
    mismatch_positions: list[int] = field(default_factory=lambda: [1])


_BACKBONES: list[Backbone] | None = None

#: Spellings accepted for the CSV's boolean columns. Anything else — including
#: a missing column — reads as False.
_TRUE_WORDS = frozenset({"true", "1", "yes"})


def _truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in _TRUE_WORDS


def _parse_mismatch_positions(value: str | None) -> list[int]:
    """Parse a backbone row's `mismatch_positions` column.

    Same comma/space-separated, 1-based convention as
    `oligo.parse_mismatch_positions` — including its standalone ``0`` means
    "no mismatch" sentinel, returning ``[]`` — but without that function's
    ``1..guide_length`` range check: a backbone isn't tied to one guide
    length, so that check happens later, against the guide actually being
    designed (`cli._resolve_mismatch_positions`, and `parse_mismatch_positions`
    itself when the TUI reads the field it seeds from this column).
    `oligo.py` already imports from this module, so it can't be imported
    back here — this is a small, deliberately narrower duplicate of just
    the two things a hand-edited CSV cell can get wrong.

    A missing or blank column (a database written before this column
    existed) defaults to ``[1]``, the mismatch position every backbone
    implicitly used before this feature existed.
    """
    stripped = (value or "").strip()
    if not stripped:
        return [1]
    try:
        positions = [int(part) for part in stripped.replace(",", " ").split()]
    except ValueError as exc:
        raise ValueError(
            f"mismatch_positions must be a comma/space-separated list of "
            f"integers; got {value!r}."
        ) from exc
    if positions == [0]:
        return []
    if 0 in positions:
        raise ValueError(
            "mismatch_positions: position 0 (no mismatch) cannot be "
            "combined with other positions."
        )
    return positions


def _bundled_data_path() -> Path:
    return Path(__file__).parent / "data" / "backbones.csv"


def _user_data_path() -> Path:
    return Path.home() / ".config" / "knockknock" / "backbones.csv"


def _bundled_backbone_names() -> set[str]:
    """Names shipped with the package — these are not the user's to delete."""
    with open(_bundled_data_path(), newline="", encoding="utf-8") as f:
        return {row["name"].strip() for row in csv.DictReader(f) if row["name"].strip()}


def bundled_backbones_path() -> Path:
    """Path to the backbones CSV shipped with the package (read-only)."""
    return _bundled_data_path()


def backbones_data_path() -> Path:
    """Path to the backbones CSV currently in effect: the user's customized
    copy if one exists, `KNOCKKNOCK_BACKBONES_PATH` if set to an existing
    file, or otherwise the bundled defaults (until the user's first save,
    which creates the customized copy at `_user_data_path()`).
    """
    return _data_path()


def _data_path() -> Path:
    """Return the user-editable backbone database when one exists.

    Installed package directories are commonly read-only. New user backbones
    therefore live in the platform-standard configuration directory, seeded
    from the bundled defaults on first save.
    """
    configured = os.environ.get("KNOCKKNOCK_BACKBONES_PATH")
    if configured:
        configured_path = Path(configured).expanduser()
        return configured_path if configured_path.exists() else _bundled_data_path()
    user_path = _user_data_path()
    return user_path if user_path.exists() else _bundled_data_path()


def load_backbones() -> list[Backbone]:
    global _BACKBONES
    if _BACKBONES is not None:
        return _BACKBONES
    backbones = []
    with open(_data_path(), newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for lineno, row in enumerate(reader, start=2):
            name = row["name"].strip()
            if not name:
                continue
            try:
                five = parse_overhang(
                    row["5_prime_site"],
                    enzyme=row["5_prime_enzyme"].strip(),
                    junction=LEFT,
                )
                three = parse_overhang(
                    row["3_prime_site"],
                    enzyme=row["3_prime_enzyme"].strip(),
                    junction=RIGHT,
                )
                mismatch_positions = _parse_mismatch_positions(
                    row.get("mismatch_positions")
                )
            except (CutSiteError, ValueError) as exc:
                raise CutSiteError(
                    f"{_data_path().name} line {lineno} ({name}): {exc}"
                ) from exc
            backbones.append(
                Backbone(
                    name=name,
                    five_prime_enzyme=row["5_prime_enzyme"].strip(),
                    five_prime_overhang=five,
                    five_prime_flank=row["5_prime_flank"].strip().upper(),
                    loop=row["loop"].strip().upper(),
                    three_prime_flank=row["3_prime_flank"].strip().upper(),
                    three_prime_enzyme=row["3_prime_enzyme"].strip(),
                    three_prime_overhang=three,
                    guide_on_5p=_truthy(row.get("guide_on_5p")),
                    # `.get`, not `[...]`: a database written before this
                    # column existed (the user's own copy under ~/.config is
                    # seeded once and then never rewritten wholesale) has no
                    # such key, and must keep loading rather than KeyError.
                    keep_sites=_truthy(row.get("keep_sites")),
                    mismatch_positions=mismatch_positions,
                )
            )
    _BACKBONES = backbones
    return backbones


def backbone_label(backbone: Backbone) -> str:
    """Human-readable backbone name with its cut enzymes appended, e.g.
    ``"pINDUCER (XhoI/EcoRI)"`` (5' enzyme, then 3').

    `Backbone.name` (the CSV/database key used for lookup) is the bare name;
    the enzymes live in their own columns so the two can't drift out of sync.
    This is the single place that recombines them for display — anywhere a
    backbone name is shown to a user should call this rather than reading
    `backbone.name` directly.
    """
    five, three = backbone.five_prime_enzyme.strip(), backbone.three_prime_enzyme.strip()
    if not five and not three:
        return backbone.name
    enzymes = five if five == three else f"{five}/{three}"
    return f"{backbone.name} ({enzymes})"


def get_backbone(name: str) -> Backbone | None:
    for b in load_backbones():
        if b.name == name:
            return b
    return None


def list_backbone_names() -> list[str]:
    return [b.name for b in load_backbones()]


_FIELDNAMES = [
    "name",
    "5_prime_enzyme",
    "5_prime_site",
    "5_prime_flank",
    "loop",
    "3_prime_flank",
    "3_prime_enzyme",
    "3_prime_site",
    "guide_on_5p",
    "keep_sites",
    "mismatch_positions",
]


def _as_row(backbone: Backbone) -> dict[str, str]:
    """`Backbone` -> one CSV row, in `_FIELDNAMES` order."""
    return {
        "name": backbone.name,
        "5_prime_enzyme": backbone.five_prime_enzyme,
        "5_prime_site": backbone.five_prime_overhang.spec,
        "5_prime_flank": backbone.five_prime_flank,
        "loop": backbone.loop,
        "3_prime_flank": backbone.three_prime_flank,
        "3_prime_enzyme": backbone.three_prime_enzyme,
        "3_prime_site": backbone.three_prime_overhang.spec,
        "guide_on_5p": "true" if backbone.guide_on_5p else "false",
        "keep_sites": "true" if backbone.keep_sites else "false",
        "mismatch_positions": ",".join(str(p) for p in backbone.mismatch_positions),
    }


def _writable_path() -> Path:
    """The database file to modify, seeding it from the bundled defaults.

    Installed package directories are commonly read-only, so the first write
    copies the shipped CSV to the user's config location and edits that. Shared
    by `save_backbone` and `delete_backbone` so they can't disagree about which
    file they are writing.
    """
    configured = os.environ.get("KNOCKKNOCK_BACKBONES_PATH")
    path = Path(configured).expanduser() if configured else _data_path()
    if path == _bundled_data_path() and not configured:
        path = _user_data_path()
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(_bundled_data_path(), path)
    return path


def save_backbone(backbone: Backbone) -> None:
    """Append `backbone` to the user-editable backbone database.

    Raises ``ValueError`` if a backbone with the same name already exists —
    callers wanting to replace one should `delete_backbone` first.
    """
    global _BACKBONES
    if any(b.name == backbone.name for b in load_backbones()):
        raise ValueError(f"A backbone named {backbone.name!r} already exists.")

    path = _writable_path()
    write_header = path.stat().st_size == 0
    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(_as_row(backbone))

    _BACKBONES = None
    load_backbones()


def delete_backbone(name: str) -> None:
    """Remove the backbone called ``name`` from the user's database.

    Rewrites the whole CSV without that row, so a mistyped custom backbone can
    be removed from the TUI instead of the user having to find and hand-edit the
    file.

    Any backbone can be deleted, bundled defaults included. The deletion is
    local: it edits the user's copy of the database, never the read-only file
    inside the installed package, so `restore_bundled_backbones` can always put
    the shipped set back.

    Raises ``ValueError`` if no backbone of that name exists.
    """
    global _BACKBONES
    if not any(b.name == name for b in load_backbones()):
        raise ValueError(f"No backbone named {name!r}.")

    path = _writable_path()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        kept = [row for row in reader if row["name"].strip() != name]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)
        writer.writeheader()
        writer.writerows(kept)

    _BACKBONES = None
    load_backbones()


def restore_bundled_backbones() -> list[str]:
    """Put back any bundled backbone the user has deleted.

    The counterpart to `delete_backbone` being allowed to remove anything: the
    shipped definitions are never actually lost, because the package's own copy
    is read-only and untouched. Backbones the user saved themselves are left
    exactly as they are.

    Returns the names that were restored, newest-last; empty if nothing was
    missing.
    """
    global _BACKBONES
    present = {b.name for b in load_backbones()}
    path = _writable_path()

    with open(_bundled_data_path(), newline="", encoding="utf-8") as f:
        missing = [
            row
            for row in csv.DictReader(f)
            if row["name"].strip() and row["name"].strip() not in present
        ]
    if not missing:
        return []

    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)
        if path.stat().st_size == 0:
            writer.writeheader()
        writer.writerows(missing)

    _BACKBONES = None
    load_backbones()
    return [row["name"].strip() for row in missing]
