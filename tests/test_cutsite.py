"""Tests for cut-site notation and the overhangs it implies."""

import pytest

from knockknock.cutsite import (
    BLUNT,
    FIVE_PRIME,
    LEFT,
    RIGHT,
    THREE_PRIME,
    CutSite,
    CutSiteError,
    Overhang,
    complement,
    ends_are_interchangeable,
    parse_cut_site,
    parse_overhang,
    reverse_complement,
)

# ---------------------------------------------------------------------------
# The symmetric-cut rule: k = L - 2c
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "spec,kind,top",
    [
        # 5' overhangs
        ("G/AATTC", FIVE_PRIME, "AATT"),  # EcoRI
        ("C/TCGAG", FIVE_PRIME, "TCGA"),  # XhoI
        ("G/GATCC", FIVE_PRIME, "GATC"),  # BamHI
        ("A/AGCTT", FIVE_PRIME, "AGCT"),  # HindIII
        ("GC/GGCCGC", FIVE_PRIME, "GGCC"),  # NotI
        ("GG/CGCC", FIVE_PRIME, "CG"),  # NarI — 2 nt
        ("CC/AGG", FIVE_PRIME, "A"),  # BstNI — 1 nt
        ("G/GACC", FIVE_PRIME, "GAC"),  # Sau96I — 3 nt, not palindromic
        ("/GATC", FIVE_PRIME, "GATC"),  # Sau3AI — cut at the 5' edge
        # blunt
        ("CCC/GGG", BLUNT, ""),  # SmaI
        ("GAT/ATC", BLUNT, ""),  # EcoRV
        ("GGA/TCC", BLUNT, ""),  # the GGN/NCC shape
        # 3' overhangs
        ("CTGCA/G", THREE_PRIME, "TGCA"),  # PstI
        ("GGTAC/C", THREE_PRIME, "GTAC"),  # KpnI
        ("GAGCT/C", THREE_PRIME, "AGCT"),  # SacI
        ("GCATG/C", THREE_PRIME, "CATG"),  # SphI
        ("CACAGT/GTG", THREE_PRIME, "AGT"),  # DraIII — 3 nt from the N region
        ("CATG/", THREE_PRIME, "CATG"),  # NlaIII — cut at the 3' edge
    ],
)
def test_overhang_derived_from_cut_position(spec, kind, top):
    """The symmetric-cut rule k = L - 2c reproduces real enzymes' overhangs.

    Covers 5', blunt and 3' cutters, edge-of-site cuts (Sau3AI, NlaIII), and
    non-palindromic asymmetric overhangs (Sau96I, DraIII) — not just the
    tidy palindromic cases.
    """
    overhang = parse_overhang(spec)
    assert overhang.kind == kind
    assert overhang.top == top


def test_same_recognition_sequence_three_different_ends():
    """CATG cut in three places gives three different ends.

    This is why the slash position, not the recognition sequence, is the
    load-bearing information: FatI, CviAII and NlaIII all read CATG.
    """
    assert parse_overhang("/CATG").describe() == "5' CATG"
    assert parse_overhang("C/ATG").describe() == "5' AT"
    assert parse_overhang("CATG/").describe() == "3' CATG"


@pytest.mark.parametrize("marker", ["/", "^", "|"])
def test_alternative_cut_markers(marker):
    """`/`, `^` and `|` are all accepted as the top-strand cut marker."""
    assert parse_overhang(f"G{marker}AATTC").top == "AATT"


def test_cut_site_geometry():
    """`CutSite` derives site/top_cut/bottom_cut/protrusion consistently for
    EcoRI, formats back to the original notation, and carries enzyme name
    and recognition sequence through to the derived `Overhang`.
    """
    site = parse_cut_site("G/AATTC", enzyme="EcoRI")
    assert site.site == "GAATTC"
    assert site.top_cut == 1
    assert site.bottom_cut == 5
    assert site.protrusion == 4
    assert site.format() == "G/AATTC"
    assert site.overhang().enzyme == "EcoRI"
    assert site.overhang().recognition == "GAATTC"


def test_palindromic_sites_give_self_complementary_overhangs():
    """A sanity property: every palindromic site's overhang pairs with itself."""
    for spec in ("G/AATTC", "C/TCGAG", "GC/GGCCGC", "CTGCA/G", "GGTAC/C", "CATG/"):
        overhang = parse_overhang(spec)
        assert reverse_complement(overhang.top) == overhang.top, spec


# ---------------------------------------------------------------------------
# Explicit overhang forms (Type IIS and other outside-cutters)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "spec,kind,top",
    [
        ("5'AATT", FIVE_PRIME, "AATT"),
        ("5'-AATT", FIVE_PRIME, "AATT"),
        ("5' AATT", FIVE_PRIME, "AATT"),
        ("3'TGCA", THREE_PRIME, "TGCA"),
        ("3'-TGCA", THREE_PRIME, "TGCA"),
        ("blunt", BLUNT, ""),
        ("BLUNT", BLUNT, ""),
        ("", BLUNT, ""),
    ],
)
def test_explicit_overhang_forms(spec, kind, top):
    """Directly-stated overhangs (for Type IIS/outside-cutters) parse across
    whitespace/hyphen spelling variants, and — unlike site-derived overhangs
    — carry no recognition sequence, since there is none to scan for.
    """
    overhang = parse_overhang(spec)
    assert overhang.kind == kind
    assert overhang.top == top
    assert overhang.recognition == ""  # nothing to scan the insert for


@pytest.mark.parametrize(
    "junction,kind,top",
    [
        # The backbone's bottom strand protrudes -> it runs 5'->3' right to
        # left, so the top-strand coordinates are the reverse complement.
        (LEFT, FIVE_PRIME, "ATAA"),
        (RIGHT, THREE_PRIME, "ATAA"),
        # The backbone's top strand protrudes -> stated form is already
        # top-strand, left to right.
        (LEFT, THREE_PRIME, "TTAT"),
        (RIGHT, FIVE_PRIME, "TTAT"),
    ],
)
def test_stated_overhang_describes_the_backbone_at_either_junction(junction, kind, top):
    """``5'TTAT``/``3'TTAT`` name the *backbone's* overhang, written 5'->3'.

    The insert then gets the complementary end, so the sequence typed into a
    site field must never turn up verbatim on an oligo. Since the backbone
    protrudes with opposite strands at the two junctions, the same text maps
    to two different top-strand sequences depending on which end it is at —
    which is the whole reason `parse_overhang` needs a junction.
    """
    overhang = parse_overhang(f"{kind[0]}'TTAT", junction=junction)
    assert overhang.kind == kind
    assert overhang.top == top
    assert overhang.backbone_top == "TTAT"  # what the user typed, unchanged
    assert overhang.describe() == f"{kind[0]}' TTAT"


def test_palindromic_overhangs_are_junction_agnostic():
    """Why the junction never mattered before: every symmetric Type II site's
    overhang is self-complementary, so left and right readings coincide.
    """
    for spec in ("5'AATT", "3'TGCA", "5'TCGA", "3'GTAC"):
        left = parse_overhang(spec, junction=LEFT)
        right = parse_overhang(spec, junction=RIGHT)
        assert left.top == right.top, spec


def test_cut_site_notation_is_top_strand_at_either_junction():
    """A slashed site names the uncut plasmid's top strand, which is already
    top-strand coordinates wherever the junction sits — so unlike a stated
    overhang it is never reverse-complemented. Sau96I's GAC is not
    self-complementary, so this would show if it were.
    """
    left = parse_overhang("G/GACC", junction=LEFT)
    right = parse_overhang("G/GACC", junction=RIGHT)
    assert left.top == right.top == "GAC"
    # Backbone-side, though, the two ends do read differently.
    assert left.backbone_top == "GTC"
    assert right.backbone_top == "GAC"


def test_interchangeable_ends_read_alike_on_the_backbone():
    """The backbone-side restatement of `ends_are_interchangeable`: a vector
    cut with one enzyme presents the same overhang at both ends, a
    directionally cut one does not.
    """
    same = (
        parse_overhang("G/AATTC", junction=LEFT),
        parse_overhang("G/AATTC", junction=RIGHT),
    )
    assert ends_are_interchangeable(*same)
    assert same[0].backbone_top == same[1].backbone_top

    directional = (
        parse_overhang("G/AATTC", junction=LEFT),
        parse_overhang("C/TCGAG", junction=RIGHT),
    )
    assert not ends_are_interchangeable(*directional)
    assert directional[0].backbone_top != directional[1].backbone_top


def test_unknown_junction_is_rejected():
    """`junction` is validated like `kind` is — a typo must not silently fall
    back to one end and design an oligo that cannot ligate.
    """
    with pytest.raises(CutSiteError, match="Unknown junction"):
        Overhang(kind=FIVE_PRIME, top="AATT", junction="middle")


def test_spec_is_retained_for_round_tripping():
    """`Overhang.spec` keeps the original text verbatim, for both site-derived
    and directly-stated forms, so it can be redisplayed or re-parsed as-is.
    """
    assert parse_overhang("G/AATTC").spec == "G/AATTC"
    assert parse_overhang("5'ACGT").spec == "5'ACGT"


# ---------------------------------------------------------------------------
# Rejections
# ---------------------------------------------------------------------------


def test_degenerate_code_inside_the_overhang_is_rejected():
    """A degenerate code is only fatal where it lands on the overhang itself.

    DraIII (``CACNNN/GTG``) is the real case: its 3' overhang *is* the N
    region, so what it leaves depends on the plasmid and no oligo can be
    ordered against it. ``C/NNGTA`` is the same shape.
    """
    with pytest.raises(CutSiteError, match="overhang itself must be concrete"):
        parse_overhang("CACNNN/GTG")
    with pytest.raises(CutSiteError, match="overhang itself must be concrete"):
        parse_overhang("C/NNGTA")


@pytest.mark.parametrize(
    "spec,top",
    [
        ("R/AATTY", "AATT"),  # ApoI — degenerate at both outer positions
        ("N/TCGAG", "TCGA"),  # an XhoI site whose first base is unrecorded
    ],
)
def test_degenerate_code_outside_the_overhang_is_accepted(spec, top):
    """Degenerate codes elsewhere in the recognition sequence are fine.

    They do not touch the bases that have to anneal, so the overhang is still
    exactly determined and the oligo is still orderable. Rejecting these was
    the old behaviour and it turned away real enzymes (ApoI) as well as
    partly-annotated plasmid maps.
    """
    overhang = parse_overhang(spec)
    assert overhang.top == top
    assert overhang.recognition == spec.replace("/", "")


def test_missing_cut_marker_is_rejected():
    """A bare recognition sequence with no `/`/`^`/`|` is ambiguous — reject
    it rather than guess where the cut falls.
    """
    with pytest.raises(CutSiteError, match="No cut position"):
        parse_overhang("GAATTC")


def test_multiple_cut_markers_are_rejected():
    """Only the top-strand cut is ever given; a second marker is nonsense
    since the bottom-strand cut is implied by the symmetric-cut rule.
    """
    with pytest.raises(CutSiteError, match="exactly one"):
        parse_overhang("G/AATT/C")


def test_unreadable_spec_is_rejected():
    """Input matching none of the accepted notations gets a clear error
    naming what forms are accepted, not a cryptic parse failure.
    """
    with pytest.raises(CutSiteError, match="Could not read"):
        parse_overhang("not a site")


def test_overhang_rejects_inconsistent_construction():
    """`Overhang.__post_init__` enforces kind/sequence consistency directly,
    not just via the parser: blunt can't carry a sequence, a 5'/3' overhang
    can't be empty, and `kind` must be one of the three known values.
    """
    with pytest.raises(CutSiteError):
        Overhang(kind=BLUNT, top="AATT")
    with pytest.raises(CutSiteError):
        Overhang(kind=FIVE_PRIME, top="")
    with pytest.raises(CutSiteError):
        Overhang(kind="sideways", top="AATT")


def test_cut_position_outside_the_site_is_rejected():
    """`top_cut` must fall within [0, len(site)] — a cut position past the
    end of the recognition sequence is meaningless, not just unusual.
    """
    with pytest.raises(CutSiteError, match="outside"):
        CutSite(site="GAATTC", top_cut=9)


# ---------------------------------------------------------------------------
# Strand helpers and end compatibility
# ---------------------------------------------------------------------------


def test_bottom_is_positional_complement_not_reverse_complement():
    """The drawn bottom strand must not be reversed — it stays column-aligned."""
    overhang = parse_overhang("G/GACC")  # 5' GAC, not self-complementary
    assert overhang.top == "GAC"
    assert overhang.bottom == "CTG"  # complement, left-to-right
    assert reverse_complement(overhang.top) == "GTC"  # what the oligo carries


def test_complement_and_reverse_complement():
    """`complement` is positional (no reversal); `reverse_complement` is the
    usual 3'->5' flip — and both tolerate the empty sequence.
    """
    assert complement("GAC") == "CTG"
    assert reverse_complement("GAC") == "GTC"
    assert complement("") == ""


def test_ends_are_interchangeable_same_enzyme():
    """The trivial case: an end is always interchangeable with itself."""
    assert ends_are_interchangeable(parse_overhang("G/AATTC"), parse_overhang("G/AATTC"))


def test_ends_are_interchangeable_compatible_enzymes():
    """BamHI and BglII leave the same GATC 5' overhang."""
    assert ends_are_interchangeable(parse_overhang("G/GATCC"), parse_overhang("A/GATCT"))


def test_blunt_ends_are_always_interchangeable():
    """Blunt ends carry no sequence to mismatch on, so any two blunt cutters
    (SmaI, EcoRV) are interchangeable — this is why blunt cloning is never
    directional.
    """
    assert ends_are_interchangeable(parse_overhang("CCC/GGG"), parse_overhang("GAT/ATC"))


def test_ends_of_different_kinds_are_not_interchangeable():
    """SacI's 3' AGCT and HindIII's 5' AGCT read alike but cannot ligate."""
    assert not ends_are_interchangeable(
        parse_overhang("A/AGCTT"), parse_overhang("GAGCT/C")
    )


def test_distinct_ends_are_not_interchangeable():
    """EcoRI's AATT and XhoI's TCGA are both 5' overhangs but don't ligate."""
    assert not ends_are_interchangeable(
        parse_overhang("G/AATTC"), parse_overhang("C/TCGAG")
    )
