"""Tests for the ASCII duplex diagram (`knockknock.duplex_view`).

These exercise the renderer directly, without Textual's async Pilot harness —
everything in `duplex_view` is a pure function of its arguments.
"""

import re
from pathlib import Path

from knockknock import duplex_view
from knockknock.cutsite import LEFT, RIGHT, parse_overhang
from knockknock.duplex_view import (
    C_ARM,
    C_FLANK,
    C_LOOP,
    C_MM,
    C_OVH,
    _split_highlighted,
    core_parts,
    render_alignment,
    render_bottom_parts,
    render_diagram,
    render_top_parts,
)
from knockknock.oligo import OligoConfig, compute_mismatched_passenger

GUIDE = "ATGCATGCATGCATGCATGCA"  # 21 nt

FLANK5 = "TCGATCGATA"  # 10 nt
LOOP = "TTCGATCG"  # 8 nt
FLANK3 = "ATCGATCGAT"  # 10 nt

#: Rich markup emitted by the renderer, e.g. ``[#ff6666]`` / ``[/#ff6666]``.
_MARKUP_RE = re.compile(r"\[/?#[0-9a-fA-F]{6}\]")


def _config(site_5: str, site_3: str, **kwargs) -> OligoConfig:
    """A minimal config with the given junctions and short, inert body parts."""
    return OligoConfig(
        five_prime_overhang=parse_overhang(site_5, junction=LEFT),
        five_prime_flank=FLANK5,
        loop=LOOP,
        three_prime_flank=FLANK3,
        three_prime_overhang=parse_overhang(site_3, junction=RIGHT),
        **kwargs,
    )


def _arms(cfg: OligoConfig) -> tuple[str, str]:
    """The (arm1, arm2) pair GUIDE produces under this config's mismatches."""
    passenger = compute_mismatched_passenger(
        GUIDE, cfg.mismatch_positions, cfg.guide_on_5p
    )
    return cfg.arms(GUIDE, passenger)


def _plain(row: str) -> str:
    """Strip Rich colour markup, leaving the drawn characters."""
    return _MARKUP_RE.sub("", row)


def _draw(cfg: OligoConfig) -> list[str]:
    """Render *cfg* against GUIDE and return the four rows as plain text."""
    arm1, arm2 = _arms(cfg)
    return [_plain(row) for row in render_diagram(cfg, arm1, arm2)]


def _colour_per_base(parts: list[tuple[str, str]]) -> list[str]:
    """Flatten colour-tagged segments into one colour entry per base."""
    return [colour for colour, seq in parts for _ in seq]


def _complement(seq: str) -> str:
    """Local positional complement, to keep the expectation independent of cutsite."""
    return "".join({"A": "T", "T": "A", "G": "C", "C": "G"}[b] for b in seq)


# ---------------------------------------------------------------------------
# Module contract
# ---------------------------------------------------------------------------


def test_duplex_view_does_not_import_textual():
    """The renderer must stay usable outside the TUI.

    That is the whole point of the module living apart from `oligo_tui`: it
    can be unit-tested without a Pilot harness and reused by a front end that
    is not a terminal.
    """
    source = Path(duplex_view.__file__).read_text(encoding="utf-8")
    assert "import textual" not in source
    assert "from textual" not in source


def test_render_diagram_returns_four_rows():
    """`render_diagram` returns exactly four strings: top, pairing, bottom, labels."""
    cfg = _config("G/AATTC", "C/TCGAG")
    rows = render_diagram(cfg, *_arms(cfg))
    assert isinstance(rows, list)
    assert len(rows) == 4
    assert all(isinstance(row, str) for row in rows)


def test_render_diagram_rows_have_no_trailing_whitespace():
    """Every row is right-stripped, so the preview has no ragged trailing blanks."""
    for row in _draw(_config("G/AATTC", "C/TCGAG")):
        assert row == row.rstrip()


# ---------------------------------------------------------------------------
# 5' overhangs — EcoRI (G/AATTC) and XhoI (C/TCGAG)
# ---------------------------------------------------------------------------


def test_five_prime_ends_put_both_overhangs_where_they_protrude():
    """With 5'-protruding ends the insert's overhangs sit on opposite strands.

    EcoRI's ``AATT`` leads the top strand at the left junction; XhoI's ``TCGA``
    protrudes on the *bottom* strand at the right junction, so it trails the
    bottom row as its complement ``AGCT``.
    """
    top, _align, bottom, _labels = _draw(_config("G/AATTC", "C/TCGAG"))
    assert "5' AATT" in top
    assert " 3'" in top
    assert top.rstrip().endswith("TCGA═════")
    assert "AGCT 5'" in bottom  # XhoI's TCGA, complemented


def test_five_prime_ends_backbone_wings_carry_the_other_strand():
    """A 5'-protruding insert leaves the backbone protruding with the *other* strand.

    So the left wing sticks out along its bottom strand (``TTAA``, the
    complement of EcoRI's ``AATT``) and the right wing along its top strand
    (XhoI's arm really is ``5'-TCGAG...-3'``).
    """
    top, _align, bottom, _labels = _draw(_config("G/AATTC", "C/TCGAG"))
    assert "═════TTAA" in bottom  # left wing, bottom strand protrudes
    assert "TCGA═════" in top  # right wing, top strand protrudes
    assert "═════    " in top  # left wing top is flush (gap, no bases)


def test_five_prime_ends_label_row_names_both_overhangs():
    """Both junctions are annotated ``5'OVH``, bracketing the core's segments."""
    *_rows, labels = _draw(_config("G/AATTC", "C/TCGAG"))
    assert labels.count("5'OVH") == 2
    assert "3'OVH" not in labels
    for name in (
        "backbone left",
        "5'flank",
        "arm1",
        "loop",
        "arm2",
        "3'flank",
        "backbone right",
    ):
        assert name in labels


# ---------------------------------------------------------------------------
# 3' overhangs — PstI (CTGCA/G) and SphI (GCATG/C)
# ---------------------------------------------------------------------------


def test_three_prime_ends_swap_which_strand_protrudes():
    """A 3'-protruding end mirrors the 5' case on every strand.

    At the left junction the insert's *bottom* strand protrudes (leading with
    ``ACGT``, the complement of PstI's ``TGCA``) while the backbone's top
    strand sticks out; at the right junction the insert's top strand runs past
    the bottom, ending on SphI's ``CATG``.
    """
    top, _align, bottom, _labels = _draw(_config("CTGCA/G", "GCATG/C"))
    assert "═════TGCA" in top  # left wing, top strand protrudes
    assert "5' TCGATCGATA" in top  # insert top starts at the core
    assert "3' ACGT" in bottom  # insert bottom leads with the overhang
    assert "CATG 3'" in top
    assert "GTAC═════" in bottom  # right wing, bottom strand protrudes


def test_three_prime_ends_label_row_names_both_overhangs():
    """Both junctions are annotated ``3'OVH`` when both ends are 3'-protruding."""
    *_rows, labels = _draw(_config("CTGCA/G", "GCATG/C"))
    assert labels.count("3'OVH") == 2
    assert "5'OVH" not in labels


def test_mixed_ends_label_each_junction_by_its_own_kind():
    """A 3'/5' pair labels the left junction ``3'OVH`` and the right ``5'OVH``."""
    *_rows, labels = _draw(_config("CTGCA/G", "C/TCGAG"))
    assert labels.count("3'OVH") == 1
    assert labels.count("5'OVH") == 1
    assert labels.index("3'OVH") < labels.index("5'OVH")


def test_short_overhang_falls_back_to_a_two_column_label():
    """Overhangs under 4 nt get the bare kind, since ``5'OVH`` needs 5 columns.

    NarI (``GG/CGCC``) leaves only 2 nt, so its label is ``5'`` rather than
    ``5'OVH``.
    """
    *_rows, labels = _draw(_config("GG/CGCC", "C/TCGAG"))
    assert labels.count("5'OVH") == 1  # only XhoI's 4-nt end is wide enough
    assert "5'" in labels


# ---------------------------------------------------------------------------
# Blunt ends
# ---------------------------------------------------------------------------


def test_blunt_ends_draw_bare_backbone_wings():
    """Blunt junctions leave nothing single-stranded, so both wings are plain ``═``."""
    top, _align, bottom, _labels = _draw(_config("blunt", "blunt"))
    assert top.strip().startswith("═════")
    assert bottom.strip().startswith("═════")
    assert top.rstrip().endswith("═════")
    assert bottom.rstrip().endswith("═════")


def test_blunt_ends_have_no_overhang_labels():
    """With no protrusion at either end there is no overhang to annotate."""
    *_rows, labels = _draw(_config("blunt", "blunt"))
    assert "5'OVH" not in labels
    assert "3'OVH" not in labels
    assert "5'flank" in labels and "3'flank" in labels


def test_blunt_ends_start_both_strands_in_the_same_column():
    """Both strands span the full width, so the two rows line up exactly."""
    cfg = _config("blunt", "blunt")
    top, align, bottom, _labels = _draw(cfg)
    assert top.index("5' ") == bottom.index("3' ")
    assert len(align.strip()) == len(cfg.core(*_arms(cfg)))


# ---------------------------------------------------------------------------
# Guide orientation
# ---------------------------------------------------------------------------


def test_guide_on_5p_false_puts_the_passenger_on_arm1():
    """By default the guide occupies arm2, so arm1 carries the mismatched passenger."""
    cfg = _config("G/AATTC", "C/TCGAG", guide_on_5p=False)
    arm1, arm2 = _arms(cfg)
    assert arm2 == GUIDE
    assert arm1 != GUIDE
    top = _draw(cfg)[0]
    assert FLANK5 + arm1 + LOOP + arm2 + FLANK3 in top


def test_guide_on_5p_true_puts_the_guide_on_arm1():
    """Flipping `guide_on_5p` swaps the two arms in the drawn top strand."""
    cfg = _config("G/AATTC", "C/TCGAG", guide_on_5p=True)
    arm1, arm2 = _arms(cfg)
    assert arm1 == GUIDE
    assert arm2 != GUIDE
    top = _draw(cfg)[0]
    assert FLANK5 + GUIDE + LOOP + arm2 + FLANK3 in top


def test_guide_orientation_changes_the_drawn_strand():
    """The two orientations render different diagrams for the same guide."""
    on_3p = _draw(_config("G/AATTC", "C/TCGAG", guide_on_5p=False))
    on_5p = _draw(_config("G/AATTC", "C/TCGAG", guide_on_5p=True))
    assert on_3p[0] != on_5p[0]
    assert on_3p[3] == on_5p[3]  # segment labels are geometry-only, so unchanged


# ---------------------------------------------------------------------------
# Mismatch highlighting
# ---------------------------------------------------------------------------


def test_mismatch_highlighted_on_both_arms_with_guide_on_arm2():
    """A mismatch marks the passenger base *and* the guide base it should pair with.

    Both halves of the bulge are coloured: passenger position 1 sits at the
    start of arm1, and the guide base opposite it is position ``len(arm1)`` on
    arm2 — i.e. the last base of the hairpin's second arm.
    """
    cfg = _config("G/AATTC", "C/TCGAG", guide_on_5p=False, mismatch_positions=[1])
    colours = _colour_per_base(core_parts(cfg, *_arms(cfg)))
    highlighted = [i for i, c in enumerate(colours) if c == C_MM]
    arm1_start = len(FLANK5)  # 10
    arm2_start = arm1_start + 21 + len(LOOP)  # 39
    assert highlighted == [arm1_start, arm2_start + 20]


def test_mismatch_highlighted_on_both_arms_with_guide_on_arm1():
    """The highlight lands at the *same columns* when the arms are swapped.

    `mismatch_positions` are counted from the base of the stem, not from
    whichever arm happens to carry the guide — arm1's front is always the
    base of the stem and arm2's back is always the base of the stem
    (arm2's front opens onto the loop instead), regardless of
    `guide_on_5p`. So flipping the guide onto arm1 changes which letters
    sit at those columns, not which columns light up.
    """
    cfg = _config("G/AATTC", "C/TCGAG", guide_on_5p=True, mismatch_positions=[1])
    colours = _colour_per_base(core_parts(cfg, *_arms(cfg)))
    highlighted = [i for i, c in enumerate(colours) if c == C_MM]
    arm1_start = len(FLANK5)  # 10
    arm2_start = arm1_start + 21 + len(LOOP)  # 39
    assert highlighted == [arm1_start, arm2_start + 20]


def test_multiple_mismatches_are_all_highlighted():
    """Every requested position is marked, on both the passenger and guide arm."""
    cfg = _config("G/AATTC", "C/TCGAG", guide_on_5p=False, mismatch_positions=[1, 5, 9])
    colours = _colour_per_base(core_parts(cfg, *_arms(cfg)))
    highlighted = [i for i, c in enumerate(colours) if c == C_MM]
    arm1_start, arm2_start = len(FLANK5), len(FLANK5) + 21 + len(LOOP)
    assert highlighted == sorted(
        [arm1_start + 0, arm1_start + 4, arm1_start + 8]
        + [arm2_start + 20, arm2_start + 16, arm2_start + 12]
    )


def test_mismatch_colour_reaches_the_rendered_row():
    """The highlight survives into the markup `render_diagram` emits."""
    cfg = _config("G/AATTC", "C/TCGAG", mismatch_positions=[1])
    top = render_diagram(cfg, *_arms(cfg))[0]
    assert f"[{C_MM}]" in top


def test_core_parts_colour_each_segment_distinctly():
    """Flanks, arms and loop each get their own colour, in 5'->3' order."""
    cfg = _config("G/AATTC", "C/TCGAG", mismatch_positions=[1])
    colours = _colour_per_base(core_parts(cfg, *_arms(cfg)))
    assert colours[0] == C_FLANK  # 5' flank
    assert colours[len(FLANK5) + 1] == C_ARM  # arm1, past the mismatch
    assert colours[len(FLANK5) + 21] == C_LOOP  # loop
    assert colours[-1] == C_FLANK  # 3' flank


def test_split_highlighted_marks_only_in_range_positions():
    """Positions outside 1..len(seq) are ignored rather than raising."""
    parts = _split_highlighted("ACGT", C_ARM, C_MM, [0, 2, 99])
    assert parts == [(C_ARM, "A"), (C_MM, "C"), (C_ARM, "GT")]


def test_split_highlighted_without_positions_is_one_segment():
    """An unmarked sequence stays a single default-coloured run."""
    assert _split_highlighted("ACGT", C_ARM, C_MM, []) == [(C_ARM, "ACGT")]
    assert _split_highlighted("", C_ARM, C_MM, [1]) == []


# ---------------------------------------------------------------------------
# Strand composition
# ---------------------------------------------------------------------------


def test_top_parts_take_the_left_5p_and_right_3p_overhangs():
    """Overhangs join the top strand only where the top strand is the one that runs on."""
    cfg = _config("G/AATTC", "GCATG/C")  # left 5' AATT, right 3' CATG
    parts = render_top_parts(cfg, core_parts(cfg, *_arms(cfg)))
    assert parts[0] == (C_OVH, "AATT")
    assert parts[-1] == (C_OVH, "CATG")


def test_bottom_parts_take_the_left_3p_and_right_5p_overhangs():
    """The mirror case: the bottom strand carries the other two protrusions.

    Bottom-strand segments are the positional complement of the top strand at
    every column, because overhangs are stored in top-strand coordinates.
    """
    cfg = _config("CTGCA/G", "C/TCGAG")  # left 3' TGCA, right 5' TCGA
    parts = render_bottom_parts(cfg, core_parts(cfg, *_arms(cfg)))
    assert parts[0] == (C_OVH, "ACGT")  # complement of TGCA
    assert parts[-1] == (C_OVH, "AGCT")  # complement of TCGA


def test_blunt_ends_add_nothing_to_either_strand():
    """Blunt junctions contribute no overhang segment to top or bottom."""
    cfg = _config("blunt", "blunt")
    core = core_parts(cfg, *_arms(cfg))
    assert render_top_parts(cfg, core) == core
    assert [seq for _c, seq in render_bottom_parts(cfg, core)] == [
        _complement(seq) for _c, seq in core
    ]


# ---------------------------------------------------------------------------
# Alignment row
# ---------------------------------------------------------------------------


def test_render_alignment_marks_paired_columns():
    """A '|' is drawn wherever the two strands are present and complementary."""
    assert render_alignment("AT", 0, "TA", 0, 2) == "||"


def test_render_alignment_skips_unpaired_columns():
    """Columns where both strands are present but do not pair stay blank."""
    assert render_alignment("AA", 0, "TG", 0, 2) == "| "


def test_render_alignment_skips_single_stranded_columns():
    """Overhang columns — only one strand present — are never marked."""
    # Bottom starts two columns in, so columns 0-1 are top-strand only.
    assert render_alignment("ATGC", 0, "CG", 2, 4) == "  ||"
    # Top starts two columns in, so columns 0-1 are bottom-strand only.
    assert render_alignment("AT", 2, "CGTA", 0, 4) == "  ||"


def test_alignment_row_covers_exactly_the_double_stranded_core():
    """Sticky ends leave the overhang columns unmarked, the core fully marked.

    The deliberate passenger mismatch does not break pairing here: it is a
    mismatch between the hairpin's two *arms*, while the bottom strand drawn
    below is the positional complement of the top strand at every column.
    """
    cfg = _config("G/AATTC", "C/TCGAG")
    arm1, arm2 = _arms(cfg)
    align = _plain(render_diagram(cfg, arm1, arm2)[1])
    assert align.count("|") == len(cfg.core(arm1, arm2))
    assert "|" not in align[: len(cfg.five_prime_overhang)]


def test_alignment_row_is_contiguous():
    """The paired region is one unbroken run, not scattered marks."""
    cfg = _config("CTGCA/G", "GCATG/C")
    align = _plain(render_diagram(cfg, *_arms(cfg))[1])
    marked = align.strip()
    assert set(marked) == {"|"}
