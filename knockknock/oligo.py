"""Oligo assembly engine — build, anneal, and batch-export oligo pairs."""

import csv
import io
import json
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from typing import Protocol

from Bio.Seq import MutableSeq, Seq
from Bio.SeqUtils import MeltingTemp as _mt

from knockknock.backbone import Backbone, backbone_label
from knockknock.cutsite import (
    FIVE_PRIME,
    THREE_PRIME,
    Overhang,
    complement,
    ends_are_interchangeable,
    invalid_dna_chars,
    reverse_complement,
    site_occurs_in,
)

# Deterministic transversion: each base maps to a non-pairing partner.
# Used to introduce the seed-region bulge on the passenger arm.
# Reproducible across runs, no RNG.
MISMATCH_TABLE: dict[str, str] = {
    "G": "T",
    "T": "G",
    "A": "C",
    "C": "A",
    "U": "G",
}


@dataclass
class OligoConfig:
    #: End the insert presents at its 5' (left) junction.
    five_prime_overhang: Overhang
    five_prime_flank: str
    loop: str
    three_prime_flank: str
    #: End the insert presents at its 3' (right) junction.
    three_prime_overhang: Overhang
    guide_on_5p: bool = False
    mismatch_positions: list[int] = field(default_factory=lambda: [1])
    backbone_name: str = "Custom"
    #: Carry the rest of each recognition sequence on the insert, so ligation
    #: rebuilds the restriction site instead of scarring it. See `site_bases`.
    keep_sites: bool = False

    @classmethod
    def from_backbone(
        cls,
        backbone: Backbone,
        mismatches: list[int] | None = None,
        keep_sites: bool | None = None,
    ) -> "OligoConfig":
        """Build a config from *backbone*'s stored design.

        ``keep_sites=None`` (the default) takes the backbone's own setting —
        it is a column of the database, since whether a vector's sites are
        normally regenerated is a property of that vector. Pass True/False to
        override it for one run, the way ``--keep-sites``/``--no-keep-sites``
        do; the distinction matters because a plain ``False`` default would
        silently switch the flag off for every backbone that stores True.

        ``mismatches=None`` (the default) takes the backbone's own
        `mismatch_positions` column the same way — a lone-``0``-disabled
        empty list included, since that's a deliberate choice and not
        "unset". Pass an explicit list to override it for one run, the way
        ``-m``/``--mismatch`` does.
        """
        return cls(
            five_prime_overhang=backbone.five_prime_overhang,
            five_prime_flank=backbone.five_prime_flank,
            loop=backbone.loop,
            three_prime_flank=backbone.three_prime_flank,
            three_prime_overhang=backbone.three_prime_overhang,
            guide_on_5p=backbone.guide_on_5p,
            mismatch_positions=(
                mismatches if mismatches is not None else backbone.mismatch_positions
            ),
            backbone_name=backbone_label(backbone),
            keep_sites=backbone.keep_sites if keep_sites is None else keep_sites,
        )

    def site_bases(self, overhang: Overhang) -> str:
        """Recognition bases this insert carries at *overhang*'s junction.

        Empty unless `keep_sites` is on and the junction actually has bases the
        insert can supply — a directly stated overhang has no recognition
        sequence to rebuild, and a degenerate one may not pin down the insert's
        own contribution. `check_sites` reports both cases rather than letting
        them pass silently.

        **Never guesses.** A degenerate code in the insert's own share makes
        this return nothing at all — it does not substitute a base, and it does
        not pad to preserve the site's length. A code on the *backbone's* share
        (``CN/TCAGC`` at a left junction, where the insert only owes ``GC``) is
        no obstacle at all, since the insert never has to synthesize it.

        These bases are fully double-stranded, so they belong to the core: they
        extend the annealed overlap and both oligos, and only the overhangs stay
        single-stranded.
        """
        if not self.keep_sites or not overhang.insert_site_bases_are_defined:
            return ""
        return overhang.insert_site_bases

    def core(self, arm1: str, arm2: str) -> str:
        """The fully double-stranded body of the insert, top strand 5'->3'.

        The recognition-sequence remainders bracket the flanks when
        `keep_sites` is on: they sit between the overhang and the flank,
        because that is where the enzyme cut them away.
        """
        return (
            self.site_bases(self.five_prime_overhang)
            + self.five_prime_flank
            + arm1
            + self.loop
            + arm2
            + self.three_prime_flank
            + self.site_bases(self.three_prime_overhang)
        )

    def arms(self, guide: str, passenger: str) -> tuple[str, str]:
        """``(arm1, arm2)`` in 5'->3' insert order, per `guide_on_5p`.

        The single place guide/passenger get ordered onto the insert's two
        arms — `assemble_oligos` and every caller that needs the same split
        (CLI warning display, the TUI preview) should call this rather than
        re-deriving the `if guide_on_5p else` swap themselves.
        """
        if self.guide_on_5p:
            return guide, passenger
        return passenger, guide


def maybe_rna(seq: str, rna: bool) -> str:
    """T->U when ``rna`` is set, else pass through unchanged.

    Public API despite being a one-liner: this is the single place the
    conversion happens, and `export.py` (candidate exports), `tui.py` and this
    module's own `oligo_rows`/`export_oligos_fasta_text` all call it instead of each
    doing their own `.replace("T", "U")`.

    Convert field-by-field, not over a whole rendered blob: several backbone
    names (e.g. ``pTRIPZ``) contain a literal ``T`` and would be corrupted by a
    blanket text replace.
    """
    if rna:
        return seq.replace("T", "U")
    return seq


@dataclass
class OligoPair:
    name: str
    forward: str
    reverse: str
    guide: str
    passenger: str
    overlap: int
    tm: float  # Tm (°C) of the inter-oligo annealed region (overlap)


def parse_mismatch_positions(raw: str, guide_length: int) -> list[int]:
    """Parse a comma/space-separated, 1-based list of mismatch positions.

    Shared strict parser for both the CLI (`knockknock oligo -m`) and any
    other caller that wants hard validation rather than best-effort parsing.
    Raises ``ValueError`` (generic wording — callers needing a flag/field
    name in the message can catch and re-raise with their own prefix) for a
    non-integer token, an empty list, a repeated position, or a position
    outside ``1..guide_length``.

    ``0`` is a standalone sentinel meaning "no mismatch": it returns ``[]``
    (a perfect-complement passenger, `compute_mismatched_passenger` applying
    no substitutions) and cannot be combined with any other position.
    """
    try:
        positions = [int(part) for part in raw.replace(",", " ").split()]
    except ValueError as exc:
        raise ValueError(
            "must be a comma- or space-separated list of positive integers."
        ) from exc
    if not positions:
        raise ValueError(
            "must contain at least one position (choose 0 to disable mismatch)."
        )
    if positions == [0]:
        return []
    if 0 in positions:
        raise ValueError("position 0 (no mismatch) cannot be combined with other positions.")
    if len(set(positions)) != len(positions):
        raise ValueError("positions must not be repeated.")
    bad = [pos for pos in positions if not 1 <= pos <= guide_length]
    if bad:
        raise ValueError(
            f"position(s) {bad} are outside the 1–{guide_length} guide range "
            "(pass 0 alone to disable mismatch introduction)."
        )
    return positions


def compute_mismatched_passenger(
    guide: str, mismatch_positions: list[int], guide_on_5p: bool = False
) -> str:
    """Reverse-complement ``guide``, then substitute each 1-based mismatch
    position via `MISMATCH_TABLE`.

    ``mismatch_positions`` are always counted from the **base of the stem**
    (the end nearest the flanking/backbone sequence, farthest from the
    loop) — never from the passenger's own 5' end as such. Which end of the
    passenger *is* the base of the stem depends on `guide_on_5p`: the
    passenger occupies `OligoConfig.core`'s ``arm1`` (front = base of stem)
    when the guide is on the 3' arm (``guide_on_5p=False``, the default),
    and ``arm2`` (back = base of stem, front = loop side) when the guide is
    on the 5' arm. Position 1 must land at the base of the stem either way,
    so the two arrangements index from opposite ends of the passenger
    string. `duplex_view.core_parts` mirrors this same base-of-stem
    convention for the diagram's mismatch highlight — keep the two in sync.

    Public API: `engine.score_one` calls it (with the default
    ``guide_on_5p=False``, since no backbone/arm order is chosen yet at
    scoring time) to fill `Candidate.sense`, and the TUI's oligo designer
    calls it (through a memoising wrapper) for its live preview.
    Deterministic by design — see `MISMATCH_TABLE`.
    """
    passenger = MutableSeq(Seq(guide).reverse_complement())
    length = len(passenger)
    for mm_pos in mismatch_positions:
        idx = (length - mm_pos) if guide_on_5p else (mm_pos - 1)
        if not (0 <= idx < length):
            continue
        nt = passenger[idx]
        if nt in MISMATCH_TABLE:
            passenger[idx] = MISMATCH_TABLE[nt]
    return str(passenger)


def _annealed_tm(overlap_seq: str) -> float:
    """Tm (°C) of the duplex formed when forward and reverse oligos anneal.

    Computed via nearest-neighbor (SantaLucia 1998 defaults) over the overlap
    region — i.e. the fully double-stranded core, excluding the overhangs at
    either end, which stay single-stranded whatever their protrusion. The
    deliberate seed-region mismatches are ignored (their effect on Tm is small
    and consistent across candidates).
    """
    if not overlap_seq:
        return 0.0
    return float(_mt.Tm_NN(overlap_seq))


def build_strands(config: OligoConfig, arm1: str, arm2: str) -> tuple[str, str]:
    """Return ``(forward, reverse)`` oligos, each written 5'->3'.

    5' overhangs extend the **5' end** of an oligo; 3' overhangs extend its
    **3' end**. Which oligo carries which follows from the junction:

    ============  ==============================  ==============================
    junction      5'-protruding                    3'-protruding
    ============  ==============================  ==============================
    left  (5')    ``forward`` 5'-prefix ``L.top``  ``reverse`` 3'-suffix ``rc(L.top)``
    right (3')    ``reverse`` 5'-prefix ``rc(R.top)``  ``forward`` 3'-suffix ``R.top``
    ============  ==============================  ==============================

    Blunt ends contribute nothing at either junction.
    """
    left = config.five_prime_overhang
    right = config.three_prime_overhang
    core = config.core(arm1, arm2)

    forward = (
        (left.top if left.kind == FIVE_PRIME else "")
        + core
        + (right.top if right.kind == THREE_PRIME else "")
    )
    reverse = (
        (reverse_complement(right.top) if right.kind == FIVE_PRIME else "")
        + reverse_complement(core)
        + (reverse_complement(left.top) if left.kind == THREE_PRIME else "")
    )
    return forward, reverse


def default_oligo_name(guide: str) -> str:
    """Fallback name for a guide with no candidate name attached (e.g. the
    `oligo` CLI subcommand's single raw-guide path). The single place that
    decides an oligo pair's name is `assemble_oligos` — this only supplies
    what it falls back to when no name is given."""
    name_parts = [guide[:8]]
    if len(guide) > 8:
        name_parts.append("...")
    return "".join(name_parts)


def assemble_oligos(
    config: OligoConfig,
    guide: str,
    passenger: str | None = None,
    name: str | None = None,
) -> OligoPair:
    if passenger is None:
        passenger = compute_mismatched_passenger(
            guide, config.mismatch_positions, config.guide_on_5p
        )

    arm1, arm2 = config.arms(guide, passenger)

    forward, reverse = build_strands(config, arm1, arm2)

    core = config.core(arm1, arm2)
    overlap = len(core)
    tm = _annealed_tm(core)

    return OligoPair(
        name=name if name is not None else default_oligo_name(guide),
        forward=forward,
        reverse=reverse,
        guide=guide,
        passenger=passenger,
        overlap=overlap,
        tm=round(tm, 1),
    )


@dataclass(frozen=True)
class InsertDuplex:
    """Column-aligned view of the annealed insert, for drawing and inspection.

    Both strands are indexed against a shared column grid ``width`` wide.
    ``top`` reads 5'->3' (it *is* the forward oligo); ``bottom`` is drawn
    left-to-right, which for the bottom strand means 3'->5' — so it is the
    reverse of the reverse oligo, and equals the positional complement of the
    top-strand sequence at every column it covers.

    The core always begins at column ``core_offset == len(left overhang)``,
    whichever way either end protrudes.
    """

    top: str
    bottom: str
    top_offset: int
    bottom_offset: int
    core_offset: int
    width: int

    @property
    def top_end(self) -> int:
        return self.top_offset + len(self.top)

    @property
    def bottom_end(self) -> int:
        return self.bottom_offset + len(self.bottom)


def make_duplex(config: OligoConfig, arm1: str, arm2: str) -> InsertDuplex:
    """Lay the annealed insert out on a column grid.

    Geometry, per junction:

    * left 5'-protruding  -> the bottom strand starts ``k`` columns in
    * left 3'-protruding  -> the top strand starts ``k`` columns in
    * right 5'-protruding -> the bottom strand runs ``k`` columns past the top
    * right 3'-protruding -> the top strand runs ``k`` columns past the bottom
    """
    left = config.five_prime_overhang
    right = config.three_prime_overhang
    core = config.core(arm1, arm2)

    top = (
        (left.top if left.kind == FIVE_PRIME else "")
        + core
        + (right.top if right.kind == THREE_PRIME else "")
    )
    bottom = (
        (left.bottom if left.kind == THREE_PRIME else "")
        + complement(core)
        + (right.bottom if right.kind == FIVE_PRIME else "")
    )
    return InsertDuplex(
        top=top,
        bottom=bottom,
        top_offset=len(left) if left.kind == THREE_PRIME else 0,
        bottom_offset=len(left) if left.kind == FIVE_PRIME else 0,
        core_offset=len(left),
        width=len(left) + len(core) + len(right),
    )


def check_ends(config: OligoConfig) -> list[str]:
    """Warn if the insert's two ends are interchangeable.

    Depends only on the backbone's junctions, not on any particular guide.
    """
    left = config.five_prime_overhang
    right = config.three_prime_overhang
    if not ends_are_interchangeable(left, right):
        return []
    detail = (
        "both ends are blunt" if left.is_blunt else f"both ends present {left.describe()}"
    )
    return [
        f"Cloning is not directional — {detail}, so the insert can ligate "
        "in either orientation."
    ]


def check_insert(config: OligoConfig, arm1: str, arm2: str) -> list[str]:
    """Note if either backbone enzyme's recognition site also occurs in the insert.

    The backbone digest happens before the insert is ligated in, so this does not
    cut the insert — it just means a copy of that restriction site ends up inside
    the cloned insert, which may matter for later steps (e.g. re-digesting the
    finished plasmid). Both strands are covered, via `cutsite.site_occurs_in`,
    which is degenerate-aware so a partly-specified site (``N/TCGAG``) still
    matches. Only the insert itself is searched — sites spanning the junctions
    depend on vector sequence this tool does not have, except when `keep_sites`
    puts the insert's half of the site back, which is deliberate and is
    reported by `check_sites` instead.
    """
    warnings: list[str] = []
    top = make_duplex(config, arm1, arm2).top

    seen: set[str] = set()
    for overhang in (config.five_prime_overhang, config.three_prime_overhang):
        site = overhang.recognition
        if not site or site in seen:
            continue
        seen.add(site)
        if site_occurs_in(site, top):
            warnings.append(
                f"The restriction sequence {site} is found in the oligo, "
                "the same restriction site used for cloning will be "
                "present in the insert as well."
            )
    return warnings


def check_sites(config: OligoConfig) -> list[str]:
    """Report where `keep_sites` could *not* deliver what it promises.

    Only problems are returned — a junction that was rebuilt cleanly says
    nothing, since the restored bases are visible in both the oligo and the
    diagram. Silence on a *failed* end would be the dangerous case, though: a
    checked box that quietly did nothing at one junction is the kind of thing
    that surfaces only after the oligos have been ordered. Nothing at all is
    reported when `keep_sites` is off.
    """
    if not config.keep_sites:
        return []

    warnings: list[str] = []
    for overhang, label in (
        (config.five_prime_overhang, "5'"),
        (config.three_prime_overhang, "3'"),
    ):
        if not overhang.recognition:
            warnings.append(
                f"Nothing to restore at the {label} junction: it was given as a "
                f"bare overhang ({overhang.describe()}), which carries no "
                "recognition sequence. Give it as a cut site (e.g. 'C/TCGAG') "
                "if the site should be rebuilt."
            )
        elif not overhang.insert_site_bases_are_defined:
            warnings.append(
                f"Cannot restore {overhang.recognition} at the {label} junction: "
                f"the base(s) the insert would have to supply "
                f"({overhang.insert_site_bases}) are degenerate. Nothing was "
                "added at that end — state those positions concretely to "
                "restore this site."
            )
        elif invalid_dna_chars(overhang.recognition)[0]:
            warnings.append(
                f"The insert carries {overhang.insert_site_bases} at the {label} "
                f"junction — the rest of {overhang.recognition}, including the "
                "degenerate position(s), comes from the vector and is not part "
                "of the oligo. Whether the site re-forms therefore depends on "
                "what the vector actually carries there."
            )
    return warnings


def validate_design(config: OligoConfig, arm1: str, arm2: str) -> list[str]:
    """All warnings for an assembled insert. Empty list means all clear."""
    return check_ends(config) + check_sites(config) + check_insert(config, arm1, arm2)


@dataclass
class BatchInput:
    """Lightweight input record for `assemble_batch`.

    `engine.Candidate` already exposes `name` and `antisense`, so it can be
    passed directly without conversion (structural-typed via `_HasGuide`).
    `sense` is accepted for forward-compatibility but is not currently honored
    (the passenger is recomputed from the configured mismatch positions).
    """

    name: str
    antisense: str
    sense: str | None = None


class _HasGuide(Protocol):
    name: str
    antisense: str


def assemble_batch(
    items: Iterable[_HasGuide],
    config: OligoConfig,
) -> list[OligoPair]:
    results = []
    for c in items:
        passenger = compute_mismatched_passenger(
            c.antisense, config.mismatch_positions, config.guide_on_5p
        )
        pair = assemble_oligos(config, c.antisense, passenger, name=c.name)
        results.append(pair)
    return results


def fasta_header(pair_name: str, strand: str, backbone_name: str) -> str:
    """Header line for one strand ("fwd"/"rev") of an assembled oligo pair.

    The single place a FASTA oligo record's header is built — `oligo.py`'s
    own `export_oligos_fasta_text` and `export.to_fasta`'s oligo section both
    call it, so the two can't drift apart. The whole name (candidate name,
    strand, and backbone — including a backbone label like
    "pINDUCER (EcoRI/XhoI)", parenthesized enzymes and all) has every space
    replaced with an underscore and every parenthesis dropped, since this
    header doubles as the oligo's name on an ordering sheet, where bare
    spaces and parentheses are not safe to assume support for.
    """
    name = f"{pair_name} {strand} {backbone_name}".replace(" ", "_")
    name = name.replace("(", "").replace(")", "")
    return f">{name}"


def export_oligos_fasta_text(
    pairs: list[OligoPair], backbone_name: str, rna: bool = False
) -> str:
    lines = []
    for pair in pairs:
        lines.append(fasta_header(pair.name, "fwd", backbone_name))
        lines.append(maybe_rna(pair.forward, rna))
        lines.append(fasta_header(pair.name, "rev", backbone_name))
        lines.append(maybe_rna(pair.reverse, rna))
        lines.append("")
    return "\n".join(lines)


def oligo_rows(
    pairs: list[OligoPair], backbone_name: str, rna: bool = False
) -> list[dict]:
    """Flat per-pair records — the single source of truth for the CSV and JSON
    renderings below, so the CLI and the TUI can never drift apart on which
    fields an oligo export carries."""
    return [
        {
            "name": pair.name,
            "backbone": backbone_name,
            "forward": maybe_rna(pair.forward, rna),
            "reverse": maybe_rna(pair.reverse, rna),
            "guide": maybe_rna(pair.guide, rna),
            "passenger": maybe_rna(pair.passenger, rna),
            "overlap": pair.overlap,
            "tm": pair.tm,
        }
        for pair in pairs
    ]


def export_oligos_csv_text(
    pairs: list[OligoPair], backbone_name: str, rna: bool = False
) -> str:
    rows = oligo_rows(pairs, backbone_name, rna=rna)
    buf = io.StringIO()
    writer = csv.DictWriter(
        buf,
        fieldnames=[
            "name",
            "backbone",
            "forward",
            "reverse",
            "guide",
            "passenger",
            "overlap",
            "tm",
        ],
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


def export_oligos_json_text(
    pairs: list[OligoPair], backbone_name: str, rna: bool = False
) -> str:
    return json.dumps(oligo_rows(pairs, backbone_name, rna=rna), indent=2)


# ---------------------------------------------------------------------------
# Batch input
# ---------------------------------------------------------------------------


class BatchInputError(ValueError):
    """Raised when a batch CSV cannot be turned into `BatchInput` records."""


#: Accepted spellings for the guide column, in priority order. First match wins.
GUIDE_COLUMNS: tuple[str, ...] = (
    "antisense",
    "guide",
    "guide_sequence",
    "Guide_Sequence",
    "Guide Sequence",
)
#: Accepted spellings for the optional passenger column.
PASSENGER_COLUMNS: tuple[str, ...] = (
    "sense",
    "passenger",
    "passenger_sequence",
    "Passenger Sequence",
)
#: Accepted spellings for the optional name column.
NAME_COLUMNS: tuple[str, ...] = ("name", "Name", "shRNA")


def _first_present(header, candidates: tuple[str, ...]) -> str | None:
    return next((c for c in candidates if c in header), None)


def read_batch_csv(
    text: str,
    *,
    normalise_guide: Callable[[str, str], str],
) -> list[BatchInput]:
    """Parse a batch CSV into `BatchInput` records.

    Column sniffing and per-row validation are domain logic, so they live here
    beside `assemble_batch` rather than inside the CLI's argument handler —
    which is also what makes them testable without invoking the CLI.

    Lines beginning with ``#`` are skipped: `export.to_csv_text` writes its
    provenance block as comments above the real header, so a candidate export
    can be fed straight back in as a batch.

    ``normalise_guide`` is injected as ``(sequence, label) -> cleaned`` so the
    caller controls how a bad guide is described; it should raise ``ValueError``
    with a user-facing message.
    """
    reader = csv.DictReader(
        line for line in io.StringIO(text) if not line.startswith("#")
    )
    rows = list(reader)
    if not rows:
        raise BatchInputError("No rows found in batch CSV.")

    header = rows[0].keys()
    guide_col = _first_present(header, GUIDE_COLUMNS)
    if guide_col is None:
        raise BatchInputError(
            "Could not find guide column in CSV. Expected one of: "
            + ", ".join(repr(c) for c in GUIDE_COLUMNS)
            + "."
        )
    passenger_col = _first_present(header, PASSENGER_COLUMNS)
    name_col = _first_present(header, NAME_COLUMNS)

    items: list[BatchInput] = []
    for i, row in enumerate(rows):
        raw_guide = (row.get(guide_col) or "").strip()
        if not raw_guide:
            continue
        # +2: one for the header row, one because humans count from 1.
        guide = normalise_guide(raw_guide, f"Guide in CSV row {i + 2}")
        name = (row.get(name_col) or "").strip() if name_col else ""
        passenger = (row.get(passenger_col) or "").strip() if passenger_col else ""
        items.append(
            BatchInput(
                name=name or f"shRNA-{i + 1}",
                antisense=guide,
                sense=passenger or None,
            )
        )

    if not items:
        raise BatchInputError("No valid guide sequences found in batch CSV.")
    return items
