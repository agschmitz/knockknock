"""Tests for oligo assembly engine and backbone loading."""

import pytest
from Bio.Seq import Seq

from knockknock.backbone import (
    Backbone,
    get_backbone,
    list_backbone_names,
    load_backbones,
)
from knockknock.cutsite import LEFT, RIGHT, complement, parse_overhang, reverse_complement
from knockknock.oligo import (
    MISMATCH_TABLE,
    BatchInput,
    BatchInputError,
    OligoConfig,
    OligoPair,
    _annealed_tm,
    assemble_batch,
    assemble_oligos,
    build_strands,
    check_ends,
    check_insert,
    check_sites,
    compute_mismatched_passenger,
    export_oligos_csv_text,
    export_oligos_fasta_text,
    export_oligos_json_text,
    fasta_header,
    make_duplex,
    maybe_rna,
    parse_mismatch_positions,
    read_batch_csv,
)

GUIDE = "ATGCATGCATGCATGCATGCA"


def _config(site_5: str, site_3: str, **kwargs) -> OligoConfig:
    """A minimal config with the given junctions and short, inert body parts.

    The junctions are parsed as LEFT/RIGHT, the way every production caller
    does — a site spec means the backbone's own end, and which strand that is
    flips between the two ends.
    """
    return OligoConfig(
        five_prime_overhang=parse_overhang(site_5, junction=LEFT),
        five_prime_flank="TCGATCGATA",
        loop="TTCGATCG",
        three_prime_flank="ATCGATCGAT",
        three_prime_overhang=parse_overhang(site_3, junction=RIGHT),
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Backbone loading
# ---------------------------------------------------------------------------


def test_load_backbones_returns_list():
    """The bundled backbones.csv loads as a non-empty list of well-formed Backbones."""
    bbs = load_backbones()
    assert isinstance(bbs, list)
    assert len(bbs) >= 1
    for bb in bbs:
        assert isinstance(bb, Backbone)
        assert bb.name
        assert bb.five_prime_enzyme
        assert bb.loop


def test_load_backbones_is_cached():
    """Repeated calls return the same list object, not a fresh re-parse."""
    bbs1 = load_backbones()
    bbs2 = load_backbones()
    assert bbs1 is bbs2


def test_get_backbone_found():
    """A known backbone name resolves to its Backbone, with junction fields intact."""
    bb = get_backbone("pINDUCER")
    assert bb is not None
    assert bb.name == "pINDUCER"
    assert bb.five_prime_enzyme == "XhoI"
    assert bb.five_prime_overhang.describe() == "5' TCGA"
    assert bb.five_prime_overhang.recognition == "CTCGAG"
    assert bb.three_prime_overhang.describe() == "5' AATT"


def test_backbone_three_prime_overhang_from_site():
    """A backbone whose junctions leave 3' protrusions parses as such.

    None of the bundled backbones (data/backbones.csv) use 3'-protruding
    enzymes, so this constructs one directly rather than via get_backbone().
    """
    bb = Backbone(
        name="Custom (PstI/KpnI)",
        five_prime_enzyme="PstI",
        five_prime_overhang=parse_overhang("CTGCA/G", enzyme="PstI"),
        five_prime_flank="TCGATCGATA",
        loop="TTCGATCG",
        three_prime_flank="ATCGATCGAT",
        three_prime_enzyme="KpnI",
        three_prime_overhang=parse_overhang("GGTAC/C", enzyme="KpnI"),
        guide_on_5p=False,
    )
    assert bb.five_prime_overhang.describe() == "3' TGCA"
    assert bb.three_prime_overhang.describe() == "3' GTAC"


def test_get_backbone_not_found():
    """An unknown backbone name returns None rather than raising."""
    assert get_backbone("nonexistent") is None


def test_list_backbone_names():
    """Returns the bundled backbones' names, including the known pINDUCER/pLKO ones."""
    names = list_backbone_names()
    assert isinstance(names, list)
    assert "pINDUCER" in names
    assert "pLKO" in names


def test_backbone_field_lengths():
    """Every bundled backbone leaves 4-base overhangs at both junctions.

    Flank and loop lengths vary per real vector (e.g. pINDUCER's 19-base
    miR30 loop vs. pLKO's 6-base one), so only the overhang length is a
    dataset-wide invariant.
    """
    for bb in load_backbones():
        assert len(bb.five_prime_overhang) == 4, f"{bb.name} 5'OVH"
        assert len(bb.three_prime_overhang) == 4, f"{bb.name} 3'OVH"


def test_plko_reproduces_the_addgene_oligo_template():
    """pLKO must emit Addgene's published pLKO.1 oligo pair, base for base.

    https://www.addgene.org/protocols/plko/ specifies:

        Forward: CCGG - 21bp sense - CTCGAG - 21bp antisense - TTTTTG
        Reverse: AATTCAAAAA - 21bp sense - CTCGAG - 21bp antisense

    where `AATTCAAAAA` is EcoRI's AATT overhang followed by rc(TTTTTG). pLKO
    is a plain stem-loop rather than a miR30 cassette, so the passenger is a
    perfect complement — hence `mismatches=[]`.
    """
    sense = "TGCCTACGTTAAGCTATACAG"
    guide = str(Seq(sense).reverse_complement())

    config = OligoConfig.from_backbone(get_backbone("pLKO"), mismatches=[])
    pair = assemble_oligos(
        config, guide, compute_mismatched_passenger(guide, [], config.guide_on_5p)
    )

    assert pair.forward == "CCGG" + sense + "CTCGAG" + guide + "TTTTTG"
    assert pair.reverse == "AATTCAAAAA" + sense + "CTCGAG" + guide


def test_phrep_reproduces_the_published_h2b_hairpin():
    """pHREP must emit the duplex published as Fig 2B of Pfisterer et al. 2026
    (PLOS One, https://doi.org/10.1371/journal.pone.0342014), the paper the
    plasmid is from (Addgene #240416).

    Esp3I (Golden Gate) leaves 4-nt 5' overhangs at both ends. The backbone's
    are `CGCT` and `TGCC` read 5'->3', so the insert presents `AGCG` on the
    forward oligo's 5' end and `GGCA` on the reverse oligo's. This is the case
    that motivated giving `parse_overhang` a junction: both ends are stated
    directly and neither is palindromic, so reading `5'CGCT` as the insert's
    end instead of the backbone's silently yields an insert that cannot
    ligate.
    """
    arm5 = "ACTTGACGAACGAGTTCATGAT"  # 5' arm, mismatched (the passenger)
    arm3 = "ATCATGAACTCGTTCGTCAAGG"  # 3' arm, the guide against chicken H2B
    loop = "GCGTGAAGCCACAGATGGC"

    backbone = get_backbone("pHREP")
    assert backbone.five_prime_overhang.backbone_top == "CGCT"
    assert backbone.three_prime_overhang.backbone_top == "TGCC"

    config = OligoConfig.from_backbone(backbone, mismatches=[1])
    pair = assemble_oligos(config, arm3, compute_mismatched_passenger(arm3, [1], False))

    # The single deliberate mismatch sits at the base of the stem, T->G.
    assert pair.passenger == arm5
    assert reverse_complement(arm5)[:-1] == arm3[:-1]
    assert reverse_complement(arm5)[-1] == "T" and arm3[-1] == "G"

    core = arm5 + loop + arm3
    assert pair.forward == "AGCG" + core
    assert pair.reverse == "GGCA" + str(Seq(core).reverse_complement())


# ---------------------------------------------------------------------------
# keep_sites: rebuilding the recognition sequence across the junction
# ---------------------------------------------------------------------------


def test_keep_sites_is_off_by_default():
    """The insert carries the overhang only unless asked otherwise, so every
    existing design is unchanged."""
    config = _config("C/TCGAG", "G/AATTC")
    assert config.keep_sites is False
    assert config.site_bases(config.five_prime_overhang) == ""
    assert config.core("AAAAA", "TTTTT").startswith(config.five_prime_flank)


def test_keep_sites_restores_the_site_across_the_ligated_junction():
    """The point of the flag: after ligation the enzyme can cut again.

    pLKO's AgeI end is the honest test — `ACCGGT` is genuinely destroyed
    without the flag. (pINDUCER's is not: its 5' flank happens to start with
    the same `G` that XhoI needs, so `CTCGAG` re-forms there by luck.)
    """
    config = _config("A/CCGGT", "G/AATTC", keep_sites=True)
    left, right = config.five_prime_overhang, config.three_prime_overhang

    # The insert supplies the far side of each cut: AgeI's trailing T at the
    # left junction, EcoRI's leading G at the right one.
    assert config.site_bases(left) == "T"
    assert config.site_bases(right) == "G"

    # Ligated top strand = backbone's head + overhang + core + backbone's tail.
    core = config.core("AAAAA", "TTTTT")
    ligated = left.site_head + left.top + core + right.top + right.site_tail
    assert "ACCGGT" in ligated
    assert "GAATTC" in ligated

    # Without the flag the insert starts at its flank instead, so the site is
    # restored only if that flank happens to open with the missing base. It
    # often does — pINDUCER's 5' flank starts with XhoI's `G`, and this file's
    # own `_config` flank starts with AgeI's `T` — so pick one that does not,
    # or the negative case passes by coincidence rather than by rule.
    without = OligoConfig(
        five_prime_overhang=parse_overhang("A/CCGGT", junction=LEFT),
        five_prime_flank="GCGATCGATA",  # opens with G, not AgeI's T
        loop="TTCGATCG",
        three_prime_flank="ATCGATCGAT",
        three_prime_overhang=parse_overhang("G/AATTC", junction=RIGHT),
    )
    scarred = left.site_head + left.top + without.core("AAAAA", "TTTTT")
    assert "ACCGGT" not in scarred


def test_keep_sites_extends_both_oligos_and_the_overlap():
    """The restored bases are double-stranded, so they join the core rather
    than the overhang: both strands grow, and so does the annealed region."""
    plain = _config("C/TCGAG", "G/AATTC")
    kept = _config("C/TCGAG", "G/AATTC", keep_sites=True)

    plain_pair = assemble_oligos(plain, GUIDE)
    kept_pair = assemble_oligos(kept, GUIDE)

    # XhoI contributes G at the left junction, EcoRI G at the right.
    assert kept_pair.forward == "TCGA" + "G" + plain.core(*plain.arms(GUIDE, plain_pair.passenger)) + "G"
    assert kept_pair.overlap == plain_pair.overlap + 2
    assert len(kept_pair.reverse) == len(plain_pair.reverse) + 2


def test_keep_sites_warns_when_there_is_no_site_to_restore():
    """A directly stated overhang carries no recognition sequence — the flag
    silently doing nothing at that end is exactly what must not happen."""
    config = _config("5'AGCG", "5'TGCC", keep_sites=True)
    messages = check_sites(config)
    assert len(messages) == 2
    assert all("Nothing to restore" in m for m in messages)
    # ...and the flag really did nothing: the core is the plain one.
    assert config.core("AAAAA", "TTTTT") == _config("5'AGCG", "5'TGCC").core(
        "AAAAA", "TTTTT"
    )


@pytest.mark.parametrize("site", ["CN/TCAGC", "N/TCGAG", "R/AATTY"])
def test_keep_sites_never_invents_a_base_for_a_degenerate_position(site):
    """A degenerate code must never reach an oligo, at either junction.

    The two escape hatches are both refusals, not guesses: if the insert's
    share of the site is degenerate it contributes *nothing* (it does not
    substitute a base, and it does not pad to keep the site's length), and if
    only the backbone's share is degenerate the insert never has to synthesize
    it in the first place. `CN/TCAGC` covers both — the insert owes the
    defined `GC` at a left junction but the degenerate `CN` at a right one.
    """
    config = _config(site, site, keep_sites=True)
    pair = assemble_oligos(config, GUIDE)

    assert not set(pair.forward + pair.reverse) - set("ACGT")

    for overhang in (config.five_prime_overhang, config.three_prime_overhang):
        added = config.site_bases(overhang)
        # Whatever was added is exactly the insert's own share, verbatim.
        assert added in ("", overhang.insert_site_bases)
        if added == "":
            continue
        assert not set(added) & set("RYSWKMBDHVN")


def test_keep_sites_handles_a_partly_specified_site_per_junction():
    """`N/TCGAG` is restorable at one junction and not the other.

    Digestion splits the site: the insert owes the trailing `G` at a left
    junction (defined, so the flag works, though whether the site re-forms
    depends on the vector's real base) but the leading `N` at a right one
    (undefined, so it cannot).
    """
    config = _config("N/TCGAG", "N/TCGAG", keep_sites=True)

    assert config.site_bases(config.five_prime_overhang) == "G"
    assert config.site_bases(config.three_prime_overhang) == ""

    left_msg, right_msg = check_sites(config)
    assert "comes from the vector and is not part of the oligo" in left_msg
    assert "Cannot restore" in right_msg
    assert "Nothing was added at that end" in right_msg


# ---------------------------------------------------------------------------
# OligoConfig
# ---------------------------------------------------------------------------


def test_config_from_backbone():
    """`OligoConfig.from_backbone` copies the backbone's junctions/loop and,
    when no `mismatches` override is given, its own `mismatch_positions`
    column (`[1]` for the bundled pINDUCER row)."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    assert config.five_prime_overhang.top == "TCGA"
    assert config.loop == "TAGTGAAGCCACAGATGTA"
    assert config.guide_on_5p is False
    assert config.backbone_name == "pINDUCER (XhoI/EcoRI)"
    assert config.mismatch_positions == [1]


def test_config_from_backbone_uses_the_backbones_own_mismatch_positions():
    """A backbone that stores something other than `[1]` is honoured, not
    silently overridden by a hardcoded default."""
    bb = get_backbone("pINDUCER")
    bb.mismatch_positions = [2, 4]
    config = OligoConfig.from_backbone(bb)
    assert config.mismatch_positions == [2, 4]


def test_config_from_backbone_custom_mismatches():
    """An explicit `mismatches` list overrides the backbone's own stored
    default."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb, mismatches=[1, 4])
    assert config.mismatch_positions == [1, 4]


def test_config_from_backbone_empty_mismatches_disables_not_defaults():
    """An explicit empty `mismatches` list (position 0, "no mismatch") must
    stay empty, not fall back to `[1]` — `from_backbone` used to test
    `mismatches or [1]`, which treated `[]` as falsy and silently
    re-enabled the default mismatch."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb, mismatches=[])
    assert config.mismatch_positions == []


# ---------------------------------------------------------------------------
# Mismatch computation
# ---------------------------------------------------------------------------


def test_mismatch_changes_passenger():
    """A single mismatch position alters that base of the reverse-complement
    passenger while keeping it a valid DNA sequence of the same length."""
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])
    # Passenger is reverse complement of guide, with pos 1 changed deterministically
    rc = str(Seq(guide).reverse_complement())
    assert len(passenger) == len(guide)
    assert passenger[0] != rc[0]
    assert set(passenger).issubset({"A", "T", "G", "C"})


def test_mismatch_position_counts_from_base_of_stem_regardless_of_arm():
    """Position 1 always lands at the base of the stem (the flank-adjacent
    end), never at the loop end — so which index of the passenger string it
    mutates flips with `guide_on_5p`.

    With the guide on arm2 (`guide_on_5p=False`, the default), the
    passenger occupies arm1 and its own front *is* the base of the stem, so
    position 1 is passenger index 0. With the guide on arm1
    (`guide_on_5p=True`), the passenger occupies arm2 — whose front opens
    onto the loop instead — so the base of the stem is the passenger's
    *last* index, and position 1 must mutate that end instead.
    """
    guide = "ATGCATGCATGCATGCATGCA"
    rc = str(Seq(guide).reverse_complement())

    passenger_3p = compute_mismatched_passenger(guide, [1], guide_on_5p=False)
    assert passenger_3p[0] != rc[0]
    assert passenger_3p[1:] == rc[1:]

    passenger_5p = compute_mismatched_passenger(guide, [1], guide_on_5p=True)
    assert passenger_5p[-1] != rc[-1]
    assert passenger_5p[:-1] == rc[:-1]


def test_mismatch_multiple_positions():
    """Several mismatch positions can be applied in one call without breaking
    the passenger's length or alphabet."""
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1, 3, 5])
    assert len(passenger) == len(guide)
    assert set(passenger).issubset({"A", "T", "G", "C"})


def test_mismatch_out_of_bounds_ignored():
    """A mismatch position past the guide's length is silently skipped, not an error."""
    guide = "ATGC"
    passenger = compute_mismatched_passenger(guide, [10])
    assert len(passenger) == 4
    assert set(passenger).issubset({"A", "T", "G", "C"})


def test_mismatch_table_is_deterministic():
    """Same input must produce same output across calls (no RNG)."""
    guide = "ATGCATGCATGCATGCATGCA"
    p1 = compute_mismatched_passenger(guide, [1, 5, 10])
    p2 = compute_mismatched_passenger(guide, [1, 5, 10])
    p3 = compute_mismatched_passenger(guide, [1, 5, 10])
    assert p1 == p2 == p3


def test_mismatch_table_correct_substitution():
    """For each base in MISMATCH_TABLE, a homopolymer guide produces the expected
    substitution at the requested passenger position."""
    # Guide of all G -> rev-comp passenger of all C; mismatch at pos 1 -> C->A
    p = compute_mismatched_passenger("G" * 10, [1])
    assert p[0] == MISMATCH_TABLE["C"]  # passenger nt was C, becomes A
    # Guide of all C -> rev-comp passenger of all G; mismatch -> G->T
    p = compute_mismatched_passenger("C" * 10, [1])
    assert p[0] == MISMATCH_TABLE["G"]
    # Guide of all A -> passenger of all T; mismatch -> T->G
    p = compute_mismatched_passenger("A" * 10, [1])
    assert p[0] == MISMATCH_TABLE["T"]
    # Guide of all T -> passenger of all A; mismatch -> A->C
    p = compute_mismatched_passenger("T" * 10, [1])
    assert p[0] == MISMATCH_TABLE["A"]


def test_mismatch_table_does_not_silently_keep_original():
    """Regression for the C->C bug: every base must map to a different base."""
    for base, replacement in MISMATCH_TABLE.items():
        assert replacement != base, f"{base} maps to itself"


def test_passenger_differs_at_mismatch_position_for_all_bases():
    """For every base at passenger position 1, the result must differ from the
    plain reverse-complement (and only at the requested position)."""
    for nt in "ACGT":
        # Guide is a 22-mer starting with this nt's complement so the passenger
        # position-1 is exactly ``nt``.
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        guide = ("ATGC" * 6)[:21] + complement[nt]
        rc = str(Seq(guide).reverse_complement())
        passenger = compute_mismatched_passenger(guide, [1])
        assert passenger[0] != rc[0], f"{nt}: passenger pos 1 unchanged"
        # Other positions must equal the plain reverse-complement.
        assert passenger[1:] == rc[1:], f"{nt}: positions >1 unexpectedly changed"


# ---------------------------------------------------------------------------
# assemble_oligos
# ---------------------------------------------------------------------------


def test_assemble_oligos_guide_on_3p():
    """End-to-end assembly for a real backbone (guide on arm2): both oligos
    carry their overhang prefix and a positive overlap/Tm come out."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])

    pair = assemble_oligos(config, guide, passenger)

    assert isinstance(pair, OligoPair)
    assert pair.forward.startswith("TCGA")
    assert pair.reverse.startswith("AATT")
    assert pair.guide == guide
    assert len(pair.forward) > len(guide)
    assert len(pair.reverse) > len(guide)
    assert pair.overlap > 0
    assert pair.tm > 0


def test_assemble_oligos_guide_on_5p():
    """Type IIS backbone (overhang stated directly) with the guide on arm1.

    No bundled backbone uses guide_on_5p, so the config is built directly
    rather than via get_backbone().
    """
    config = _config("5'ACGT", "5'AAAA", guide_on_5p=True)
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])

    pair = assemble_oligos(config, guide, passenger)

    # guide should be arm1 when guide_on_5p
    assert pair.forward.startswith("ACGT")
    assert pair.forward.find(guide) < pair.forward.find(passenger)
    # guide should be in the arm1 position (right after 5'flank)
    arm1_start = len(config.five_prime_overhang) + len(config.five_prime_flank)
    fwd_arm1 = pair.forward[arm1_start : arm1_start + len(guide)]
    assert fwd_arm1 == guide


def test_assemble_oligos_auto_passenger():
    """Omitting the passenger makes `assemble_oligos` derive one itself
    (via `compute_mismatched_passenger` and the config's mismatch positions)."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"

    pair = assemble_oligos(config, guide)
    assert pair.passenger
    assert len(pair.passenger) == len(guide)


def test_assemble_oligos_overlap_length():
    """`pair.overlap` equals the full double-stranded core length (flanks +
    guide + loop + passenger), matching `OligoConfig.core`."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])

    pair = assemble_oligos(config, guide, passenger)

    # Compared against `core()` itself rather than a hand-summed segment list:
    # the core grew a `keep_sites` member at each end, and a re-derived sum
    # silently goes stale every time that happens.
    arm1, arm2 = config.arms(guide, passenger)
    assert pair.overlap == len(config.core(arm1, arm2))

    # Spelled out once, so the identity above still documents what the core is.
    assert pair.overlap == (
        len(config.site_bases(config.five_prime_overhang))
        + len(config.five_prime_flank)
        + len(guide)
        + len(config.loop)
        + len(passenger)
        + len(config.three_prime_flank)
        + len(config.site_bases(config.three_prime_overhang))
    )


# ---------------------------------------------------------------------------
# assemble_batch
# ---------------------------------------------------------------------------


def test_assemble_batch():
    """`assemble_batch` runs each `BatchInput` through `assemble_oligos`,
    preserving input order and each item's name on the resulting pair."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)

    candidates = [
        BatchInput("shRNA-1", "ATGCATGCATGCATGCATGCA"),
        BatchInput("shRNA-2", "GCTAGCTAGCTAGCTAGCTAG"),
        BatchInput("shRNA-3", "TTGATTGATTGATTGATTGAT"),
    ]

    pairs = assemble_batch(candidates, config)
    assert len(pairs) == 3
    assert pairs[0].name == "shRNA-1"
    assert pairs[1].name == "shRNA-2"
    assert pairs[2].name == "shRNA-3"
    assert all(p.forward.startswith("TCGA") for p in pairs)
    assert all(p.reverse.startswith("AATT") for p in pairs)


# ---------------------------------------------------------------------------
# fasta_header
# ---------------------------------------------------------------------------


def test_fasta_header_uses_short_strand_labels():
    """"fwd"/"rev", not "Forward"/"Reverse" — short enough for an ordering
    sheet's name column."""
    assert fasta_header("sh-1", "fwd", "Custom") == ">sh-1_fwd_Custom"
    assert fasta_header("sh-1", "rev", "Custom") == ">sh-1_rev_Custom"


def test_fasta_header_strips_spaces_and_parens_from_backbone_name():
    """A backbone label like "pINDUCER (EcoRI/XhoI)" carries a space and
    parenthesized enzymes; neither may survive into the header, since the
    whole thing doubles as the oligo's name on an ordering sheet, which is
    not safe to assume supports spaces or parentheses."""
    header = fasta_header("sh-1", "fwd", "pINDUCER (EcoRI/XhoI)")
    assert header == ">sh-1_fwd_pINDUCER_EcoRI/XhoI"
    assert " " not in header
    assert "(" not in header and ")" not in header


# ---------------------------------------------------------------------------
# export_oligos_fasta_text
# ---------------------------------------------------------------------------


def test_export_oligos_fasta_text_format():
    """A single pair renders as a 4-line FASTA-style block: fwd header +
    sequence, then rev header + sequence — the header is one space-free
    token (name_strand_backbone), including the backbone's own spaces."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    pair = assemble_oligos(config, guide)

    text = export_oligos_fasta_text([pair], "pINDUCER (EcoRI/XhoI)")
    lines = text.strip().split("\n")

    # Format: >name_fwd_backbone\nseq\n>name_rev_backbone\nseq
    assert lines[0].startswith(">")
    assert "_fwd_" in lines[0]
    assert "pINDUCER" in lines[0]
    assert " " not in lines[0] and "(" not in lines[0] and ")" not in lines[0]
    assert len(lines[1]) == len(pair.forward)
    assert lines[2].startswith(">")
    assert "_rev_" in lines[2]
    assert " " not in lines[2] and "(" not in lines[2] and ")" not in lines[2]
    assert len(lines[3]) == len(pair.reverse)


def test_export_oligos_fasta_text_multiple():
    """Each pair in a multi-pair export gets its own fwd/rev block, keyed by
    its own name."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)

    candidates = [
        BatchInput("sh-1", "ATGCATGCATGCATGCATGCA"),
        BatchInput("sh-2", "GCTAGCTAGCTAGCTAGCTAG"),
    ]
    pairs = assemble_batch(candidates, config)
    text = export_oligos_fasta_text(pairs, "Test")
    # Should have 4 FASTA entries + 1 blank line between pairs = 9 lines total (including trailing blank)
    assert text.count(">sh-1_fwd") == 1
    assert text.count(">sh-2_rev") == 1


def test_export_oligos_csv_text_content():
    """CSV export has the full fixed column set and each row's values match
    the source `OligoPair`."""
    import csv
    import io

    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    pair = assemble_oligos(config, guide)

    text = export_oligos_csv_text([pair], "pINDUCER (EcoRI/XhoI)")
    rows = list(csv.DictReader(io.StringIO(text)))

    assert rows[0].keys() == {
        "name",
        "backbone",
        "forward",
        "reverse",
        "guide",
        "passenger",
        "overlap",
        "tm",
    }
    assert len(rows) == 1
    assert rows[0]["backbone"] == "pINDUCER (EcoRI/XhoI)"
    assert rows[0]["guide"] == pair.guide
    assert rows[0]["forward"] == pair.forward
    assert rows[0]["reverse"] == pair.reverse
    assert rows[0]["overlap"] == str(pair.overlap)


def test_export_oligos_json_text_content():
    """JSON export is a list of one object per pair, with fields matching the
    source `OligoPair` (including numeric `overlap`/`tm`, not stringified)."""
    import json

    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    pair = assemble_oligos(config, guide)

    payload = json.loads(export_oligos_json_text([pair], "pINDUCER (EcoRI/XhoI)"))

    assert len(payload) == 1
    row = payload[0]
    assert row["backbone"] == "pINDUCER (EcoRI/XhoI)"
    assert row["guide"] == pair.guide
    assert row["forward"] == pair.forward
    assert row["reverse"] == pair.reverse
    assert row["overlap"] == pair.overlap
    assert row["tm"] == pair.tm


def test_export_oligos_rna_converts_per_field_not_whole_blob():
    """Regression: a whole-text T->U replace corrupts backbone names that
    contain a literal T (e.g. pTRIPZ). RNA conversion must be per-field."""
    bb = get_backbone("pTRIPZ")
    config = OligoConfig.from_backbone(bb)
    pair = assemble_oligos(config, "ATGCATGCATGCATGCATGCA")

    text = export_oligos_fasta_text([pair], "pTRIPZ (EcoRI/XhoI)", rna=True)
    assert "pTRIPZ" in text  # backbone name in the header is untouched
    assert "U" in pair_body_lines(text)[0]  # but the sequence itself is RNA
    assert "T" not in pair_body_lines(text)[0]

    rows = export_oligos_csv_text([pair], "pTRIPZ (EcoRI/XhoI)", rna=True)
    assert "pTRIPZ" in rows

    import json

    payload = json.loads(export_oligos_json_text([pair], "pTRIPZ (EcoRI/XhoI)", rna=True))
    assert payload[0]["backbone"] == "pTRIPZ (EcoRI/XhoI)"
    assert "U" in payload[0]["forward"]


def pair_body_lines(text: str) -> list[str]:
    """Non-header lines of `export_oligos_fasta_text` output, for RNA assertions."""
    return [line for line in text.splitlines() if line and not line.startswith(">")]


def test_maybe_rna():
    """``rna=False`` passes the sequence through unchanged; ``rna=True`` converts T->U."""
    assert maybe_rna("ATGCT", rna=False) == "ATGCT"
    assert maybe_rna("ATGCT", rna=True) == "AUGCU"


def test_oligo_config_arms_respects_guide_on_5p():
    """`OligoConfig.arms` orders (guide, passenger) onto (arm1, arm2) per
    `guide_on_5p`, the single ordering helper all callers share."""
    config = _config("G/AATTC", "C/TCGAG", guide_on_5p=True)
    assert config.arms("GUIDE", "PASSENGER") == ("GUIDE", "PASSENGER")

    config = _config("G/AATTC", "C/TCGAG", guide_on_5p=False)
    assert config.arms("GUIDE", "PASSENGER") == ("PASSENGER", "GUIDE")


# ---------------------------------------------------------------------------
# parse_mismatch_positions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw", ["abc", "", "1,1", "1,99", "0,1", "1,0"])
def test_parse_mismatch_positions_rejects_invalid(raw):
    """Non-numeric input, an empty list, a repeated position, an
    out-of-range position, and 0 combined with another position each raise
    ValueError."""
    with pytest.raises(ValueError):
        parse_mismatch_positions(raw, guide_length=22)


def test_parse_mismatch_positions_accepts_comma_or_space():
    """Comma- and space-separated position lists parse to the same result."""
    assert parse_mismatch_positions("1,2", 22) == [1, 2]
    assert parse_mismatch_positions("1 2", 22) == [1, 2]


def test_parse_mismatch_positions_zero_disables():
    """A standalone 0 disables mismatch introduction: it parses to an empty
    list rather than being rejected as out of the 1..guide_length range."""
    assert parse_mismatch_positions("0", 22) == []


# ---------------------------------------------------------------------------
# Overhang kinds: where the extra bases land on each oligo
# ---------------------------------------------------------------------------


def test_five_prime_overhangs_extend_the_5_prime_ends():
    """EcoRI/XhoI: both overhangs prepend, as they always have."""
    config = _config("G/AATTC", "C/TCGAG")
    core = config.core("AAAAA", "TTTTT")
    forward, reverse = build_strands(config, "AAAAA", "TTTTT")
    assert forward == "AATT" + core
    assert reverse == "TCGA" + str(Seq(core).reverse_complement())


def test_three_prime_overhangs_extend_the_3_prime_ends():
    """PstI/KpnI: 3' protrusions append instead of prepend."""
    config = _config("CTGCA/G", "GGTAC/C")
    core = config.core("AAAAA", "TTTTT")
    forward, reverse = build_strands(config, "AAAAA", "TTTTT")
    # Right end protrudes on the top strand; left end on the bottom strand.
    assert forward == core + "GTAC"
    assert reverse == str(Seq(core).reverse_complement()) + "TGCA"


def test_mixed_overhang_kinds():
    """EcoRI (5') on the left, KpnI (3') on the right — both land on forward."""
    config = _config("G/AATTC", "GGTAC/C")
    core = config.core("AAAAA", "TTTTT")
    forward, reverse = build_strands(config, "AAAAA", "TTTTT")
    assert forward == "AATT" + core + "GTAC"
    assert reverse == str(Seq(core).reverse_complement())


def test_blunt_ends_add_nothing():
    """Blunt junctions contribute no extra bases at either end: both oligos
    are exactly the core (forward) or its reverse complement (reverse)."""
    config = _config("CCC/GGG", "CCC/GGG")
    core = config.core("AAAAA", "TTTTT")
    forward, reverse = build_strands(config, "AAAAA", "TTTTT")
    assert forward == core
    assert reverse == str(Seq(core).reverse_complement())


@pytest.mark.parametrize(
    "site_5,site_3,expect",
    [
        # (forward, reverse) with the non-palindromic TTAT stated at one end.
        # Whichever end and whichever protrusion, the oligo carries ATAA —
        # rc("TTAT") — because TTAT is what the *backbone* presents.
        ("5'TTAT", "blunt", ("ATAA" + "{core}", "{rc}")),
        ("3'TTAT", "blunt", ("{core}", "{rc}" + "ATAA")),
        ("blunt", "5'TTAT", ("{core}", "ATAA" + "{rc}")),
        ("blunt", "3'TTAT", ("{core}" + "ATAA", "{rc}")),
    ],
)
def test_stated_overhang_never_lands_verbatim_on_an_oligo(site_5, site_3, expect):
    """A site field states the backbone's overhang, so the oligo gets its
    reverse complement — at both junctions and for both protrusion kinds.

    The bug this pins: ``5'TTAT`` used to put TTAT straight onto the forward
    oligo's 5' end (i.e. read as the *insert's* end) while ``3'TTAT`` put it
    on the backbone, so the same notation meant two different things. Every
    real Type II overhang is palindromic, which is why only a stated,
    non-palindromic end like TTAT can tell the two readings apart.
    """
    config = _config(site_5, site_3)
    core = config.core("AAAAA", "TTTTT")
    rc = str(Seq(core).reverse_complement())
    forward, reverse = build_strands(config, "AAAAA", "TTTTT")

    assert forward == expect[0].format(core=core, rc=rc)
    assert reverse == expect[1].format(core=core, rc=rc)


def test_stated_overhang_anneals_to_the_backbone_end_it_names():
    """The ligation check, done base by base on the duplex.

    Wherever the backbone's overhang sits, the strand the *insert* brings to
    those columns must be its positional complement — otherwise the ends do
    not anneal, which is the practical consequence of getting the reading
    direction wrong.
    """
    config = _config("5'TTAT", "3'TTAT")
    duplex = make_duplex(config, "AAAAA", "TTTTT")

    left, right = config.five_prime_overhang, config.three_prime_overhang
    assert left.backbone_top == right.backbone_top == "TTAT"

    # Left end: the backbone holds the bottom strand, the insert the top.
    assert duplex.top[: len(left)] == complement(left.bottom)
    # Right end: same junction geometry mirrored — the insert's top strand
    # runs past the backbone's protruding bottom strand.
    assert duplex.top[-len(right) :] == complement(right.bottom)
    # And read 5'->3' along the backbone's own strand, both are what was typed.
    assert reverse_complement(duplex.top[: len(left)]) == "TTAT"
    assert reverse_complement(duplex.top[-len(right) :]) == "TTAT"


def test_non_palindromic_right_overhang_is_reverse_complemented():
    """A 5'-protruding right end puts rc(top) on the reverse oligo, not top.

    Palindromic sites hide this (rc(TCGA) == TCGA), so use Sau96I's GGACC,
    whose 5' overhang GAC is not self-complementary.
    """
    config = _config("CCC/GGG", "G/GACC")
    assert config.three_prime_overhang.describe() == "5' GAC"
    _, reverse = build_strands(config, "AAAAA", "TTTTT")
    assert reverse.startswith("GTC")  # rc("GAC"), not "GAC"


# ---------------------------------------------------------------------------
# make_duplex
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "site_5,site_3",
    [
        ("G/AATTC", "C/TCGAG"),  # 5' / 5'
        ("CTGCA/G", "GGTAC/C"),  # 3' / 3'
        ("G/AATTC", "GGTAC/C"),  # 5' / 3'
        ("CTGCA/G", "C/TCGAG"),  # 3' / 5'
        ("CCC/GGG", "CCC/GGG"),  # blunt / blunt
        ("GG/CGCC", "CATG/"),  # 2 nt 5' / 4 nt 3'
    ],
)
def test_duplex_strands_match_the_oligos(site_5, site_3):
    """The drawn duplex and the assembled oligos must be the same molecule."""
    config = _config(site_5, site_3)
    passenger = compute_mismatched_passenger(GUIDE, [1])
    forward, reverse = build_strands(config, passenger, GUIDE)
    duplex = make_duplex(config, passenger, GUIDE)

    # Top row is the forward oligo; the bottom row is drawn left-to-right,
    # i.e. 3'->5', so reversing it recovers the reverse oligo.
    assert duplex.top == forward
    assert duplex.bottom[::-1] == reverse


@pytest.mark.parametrize(
    "site_5,site_3",
    [
        ("G/AATTC", "C/TCGAG"),
        ("CTGCA/G", "GGTAC/C"),
        ("G/AATTC", "GGTAC/C"),
        ("CTGCA/G", "C/TCGAG"),
        ("CCC/GGG", "CCC/GGG"),
        ("GG/CGCC", "CATG/"),
    ],
)
def test_duplex_geometry(site_5, site_3):
    """The core starts at len(left overhang) whichever way the ends protrude."""
    config = _config(site_5, site_3)
    duplex = make_duplex(config, "AAAAA", "TTTTT")
    left, right = config.five_prime_overhang, config.three_prime_overhang
    core = config.core("AAAAA", "TTTTT")

    assert duplex.core_offset == len(left)
    assert duplex.width == len(left) + len(core) + len(right)
    assert duplex.width == max(duplex.top_end, duplex.bottom_end)
    # Exactly one strand is indented on the left, and only when it protrudes.
    assert min(duplex.top_offset, duplex.bottom_offset) == 0
    assert max(duplex.top_offset, duplex.bottom_offset) == len(left)
    # Every column the two strands share is a base pair.
    overlap_start = max(duplex.top_offset, duplex.bottom_offset)
    overlap_end = min(duplex.top_end, duplex.bottom_end)
    top = duplex.top[overlap_start - duplex.top_offset : overlap_end - duplex.top_offset]
    bottom = duplex.bottom[
        overlap_start - duplex.bottom_offset : overlap_end - duplex.bottom_offset
    ]
    assert top == str(Seq(bottom).complement())
    assert overlap_end - overlap_start == len(core)


def test_duplex_left_five_prime_indents_the_bottom_strand():
    """A 5'-protruding left junction: the top strand starts at column 0 and
    the bottom strand is indented by the overhang length."""
    duplex = make_duplex(_config("G/AATTC", "CCC/GGG"), "AAAAA", "TTTTT")
    assert duplex.top_offset == 0
    assert duplex.bottom_offset == 4


def test_duplex_left_three_prime_indents_the_top_strand():
    """A 3'-protruding left junction: the bottom strand starts at column 0 and
    the top strand is indented by the overhang length (mirror of the 5' case)."""
    duplex = make_duplex(_config("CTGCA/G", "CCC/GGG"), "AAAAA", "TTTTT")
    assert duplex.top_offset == 4
    assert duplex.bottom_offset == 0


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_check_ends_flags_identical_overhangs():
    """Both junctions using the same enzyme (same overhang) triggers a
    non-directional-cloning warning."""
    warnings = check_ends(_config("G/AATTC", "G/AATTC"))
    assert len(warnings) == 1
    assert "not directional" in warnings[0]


def test_check_ends_flags_blunt_blunt():
    """Blunt ends at both junctions are also non-directional, with wording
    that specifically calls out "blunt" rather than a shared overhang."""
    warnings = check_ends(_config("CCC/GGG", "CCC/GGG"))
    assert len(warnings) == 1
    assert "blunt" in warnings[0]


def test_check_ends_flags_compatible_ends_from_different_enzymes():
    """BamHI and BglII leave the same GATC overhang, so the insert can flip."""
    assert check_ends(_config("G/GATCC", "A/GATCT"))


def test_check_ends_silent_for_directional_cloning():
    """Distinct, non-interchangeable junctions (5'/5' and 5'/3' pairs) produce no warnings."""
    assert check_ends(_config("G/AATTC", "C/TCGAG")) == []
    assert check_ends(_config("G/AATTC", "GGTAC/C")) == []


def test_check_insert_finds_internal_site():
    """A cloning enzyme's recognition site planted inside the loop is caught
    and named in the warning."""
    config = OligoConfig.from_backbone(get_backbone("pINDUCER"))
    config.loop = "GAATTCGA"  # an EcoRI site planted in the loop
    warnings = check_insert(config, "AAAAA", "TTTTT")
    assert len(warnings) == 1
    assert "The restriction sequence GAATTC is found in the oligo" in warnings[0]


def test_check_insert_finds_site_on_the_bottom_strand():
    """A non-palindromic site is caught in either orientation."""
    config = _config("CCC/GGG", "G/GACC")
    config.loop = "GGTCCAAA"  # rc("GGACC") == "GGTCC"
    warnings = check_insert(config, "AAAAA", "TTTTT")
    assert any("GGACC" in w for w in warnings)


def test_check_insert_silent_when_clean():
    """No warnings when neither enzyme's recognition site appears in the insert."""
    config = _config("G/AATTC", "C/TCGAG")
    assert check_insert(config, "AAAAA", "TTTTT") == []


def test_check_insert_silent_without_a_recognition_sequence():
    """Type IIS junctions state the overhang directly, so there is nothing to scan."""
    config = _config("5'ACGT", "3'TGCA")
    config.loop = "GAATTCGA"
    assert check_insert(config, "AAAAA", "TTTTT") == []


# ---------------------------------------------------------------------------
# Tm of the inter-oligo annealed region
# ---------------------------------------------------------------------------


def test_annealed_tm_zero_for_empty_overlap():
    """An empty overlap sequence short-circuits to 0.0 rather than calling
    into the nearest-neighbor Tm calculation."""
    assert _annealed_tm("") == 0.0


def test_annealed_tm_in_range():
    """For a representative pINDUCER-backbone overlap, Tm should fall in a
    biologically reasonable range."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])
    pair = assemble_oligos(config, guide, passenger)
    assert 30 < pair.tm < 95


# ---------------------------------------------------------------------------
# Round-trip complement
# ---------------------------------------------------------------------------


def test_forward_reverse_complement_consistency():
    """After stripping each oligo's 5' overhang, the forward and reverse
    oligos are reverse complements of each other over the annealed region —
    they describe the same double-stranded molecule."""
    bb = get_backbone("pINDUCER")
    config = OligoConfig.from_backbone(bb)
    guide = "ATGCATGCATGCATGCATGCA"
    passenger = compute_mismatched_passenger(guide, [1])

    pair = assemble_oligos(config, guide, passenger)

    # Both ends of this backbone are 5'-protruding, so each oligo carries its
    # overhang as a 5' prefix and the remainder is the shared duplex core.
    fwd_overlap = pair.forward[len(config.five_prime_overhang) :]  # skip 5'OVH
    rev_overlap = pair.reverse[len(config.three_prime_overhang) :]  # skip 3'OVH

    # reverse complement of the forward overlap should equal the reverse overlap
    fwd_rc = str(Seq(fwd_overlap).reverse_complement())
    assert fwd_rc == rev_overlap


# ---------------------------------------------------------------------------
# read_batch_csv — column sniffing, moved out of the CLI so it is testable
# without invoking Typer
# ---------------------------------------------------------------------------


def _identity_guide(seq: str, label: str) -> str:
    return seq.upper()


@pytest.mark.parametrize(
    "column", ["antisense", "guide", "guide_sequence", "Guide_Sequence", "Guide Sequence"]
)
def test_read_batch_csv_accepts_every_guide_column_spelling(column):
    """Each documented spelling of the guide column is recognised."""
    text = f"{column}\nATGCATGCATGCATGCATGCAT\n"
    items = read_batch_csv(text, normalise_guide=_identity_guide)
    assert [i.antisense for i in items] == ["ATGCATGCATGCATGCATGCAT"]


def test_read_batch_csv_skips_provenance_comment_lines():
    """`knockknock run -f csv` writes '#' metadata above the header, so its own
    output can be fed straight back in as a batch."""
    text = (
        "# knockknock_version: 0.1.0\n# gc_min: 20.0\nantisense\nATGCATGCATGCATGCATGCAT\n"
    )
    items = read_batch_csv(text, normalise_guide=_identity_guide)
    assert len(items) == 1


def test_read_batch_csv_uses_the_name_column_when_present():
    """A name column labels each pair in the export."""
    text = "name,antisense\nmy-shRNA,ATGCATGCATGCATGCATGCAT\n"
    items = read_batch_csv(text, normalise_guide=_identity_guide)
    assert items[0].name == "my-shRNA"


def test_read_batch_csv_generates_names_when_the_column_is_absent():
    """Without a name column, rows are numbered from 1."""
    text = "antisense\nATGCATGCATGCATGCATGCAT\nGCTAGCTAGCTAGCTAGCTAGC\n"
    items = read_batch_csv(text, normalise_guide=_identity_guide)
    assert [i.name for i in items] == ["shRNA-1", "shRNA-2"]


def test_read_batch_csv_skips_blank_guide_rows():
    """A trailing or empty row must not become an empty oligo pair."""
    text = "antisense\nATGCATGCATGCATGCATGCAT\n\n   \n"
    items = read_batch_csv(text, normalise_guide=_identity_guide)
    assert len(items) == 1


def test_read_batch_csv_rejects_a_file_with_no_guide_column():
    """A CSV without any recognised guide column is an error, not an empty run."""
    with pytest.raises(BatchInputError, match="Could not find guide column"):
        read_batch_csv("name,score\nfoo,1\n", normalise_guide=_identity_guide)


def test_read_batch_csv_rejects_an_empty_file():
    """No rows at all is an error."""
    with pytest.raises(BatchInputError, match="No rows found"):
        read_batch_csv("", normalise_guide=_identity_guide)


def test_read_batch_csv_rejects_a_file_whose_guide_rows_are_all_blank():
    """Rows that exist but carry no guide must not silently produce zero pairs.

    (Wholly empty lines are dropped by csv.DictReader before we see them, which
    is the "No rows found" case above — this is the distinct one where the rows
    are real but the guide cell is empty.)
    """
    with pytest.raises(BatchInputError, match="No valid guide sequences"):
        read_batch_csv(
            "antisense,name\n,first\n,second\n", normalise_guide=_identity_guide
        )


def test_read_batch_csv_reports_the_spreadsheet_row_number_on_a_bad_guide():
    """The row number in the message must match what the user sees in Excel:
    the header is row 1, so the first data row is row 2."""

    def _reject(seq, label):
        raise ValueError(f"{label} is bad.")

    with pytest.raises(ValueError, match="row 2"):
        read_batch_csv("antisense\nBOGUS\n", normalise_guide=_reject)
