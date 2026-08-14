"""Tests for core engine functions."""

import pytest
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from knockknock.engine import (
    MOTIFS_TO_AVOID,
    AnalysisParameterError,
    AnalysisParams,
    make_shrnas,
    process_shrnas,
    score_one,
    validate_analysis_params,
    validate_top_n,
)
from knockknock.reference import (
    ReferenceData,
    load_coeffs,
    load_deltaG,
    load_dG_bins,
    load_reference_data,
)


@pytest.fixture
def ref() -> ReferenceData:
    return load_reference_data()


def test_make_shrnas_basic():
    """Slides a 22-nt window across the reverse complement, one candidate per
    position, and the first/last windows land where expected."""
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    params = AnalysisParams(length=22)
    result = make_shrnas(mrnas, params)
    assert len(result["test"]) == 4
    assert result["test"][0][1] == "ATCAGCTAGCTAGCTAGCTACG"
    assert result["test"][3][1] == "AGCTAGCTAGCTAGCTACGCAT"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"gc_min": 80.0, "gc_max": 20.0}, "GC limits"),
        ({"gc_min": float("nan")}, "finite"),
        ({"gc_max": float("inf")}, "finite"),
        ({"min_efficacy_score": float("nan")}, "Minimum efficacy score"),
        ({"multitarget": 0}, "multitarget"),
        ({"mismatch": [23]}, "mismatch"),
    ],
)
def test_validate_analysis_params_rejects_invalid_shared_settings(kwargs, message):
    """Every frontend receives the same filter and analysis bounds here."""
    with pytest.raises(AnalysisParameterError, match=message):
        validate_analysis_params(**kwargs)


def test_validate_analysis_params_preserves_every_valid_setting():
    """The shared boundary returns an engine-ready copy of every option."""
    params = validate_analysis_params(
        length=20,
        gc_min=0.0,
        gc_max=100.0,
        min_efficacy_score=-1.5,
        sort_by_dg=True,
        remove_motifs=True,
        multitarget=2,
        rna_output=True,
        mismatch=[20],
        seed=42,
    )

    assert params == AnalysisParams(
        length=20,
        gc_min=0.0,
        gc_max=100.0,
        min_efficacy_score=-1.5,
        sort_by_dg=True,
        remove_motifs=True,
        multitarget=2,
        rna_output=True,
        stem_mismatch=[20],
        seed=42,
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length": 19},
        {"length": 23},
        {"gc_min": -0.1},
        {"gc_max": 100.1},
        {"gc_min": 50.0, "gc_max": 50.0},
        {"min_efficacy_score": float("inf")},
        {"min_efficacy_score": float("-inf")},
        {"mismatch": [0]},
    ],
)
def test_validate_analysis_params_rejects_all_boundaries(kwargs):
    """Bounds are exclusive only where the command help promises they are.

    A literal ``0`` inside the list (as opposed to an explicit *empty*
    list) is still just an out-of-range position here — the "0 means
    disable" sentinel is a string-level convention handled by
    `oligo.parse_mismatch_positions` before the list ever reaches this
    function.
    """
    with pytest.raises(AnalysisParameterError):
        validate_analysis_params(**kwargs)


def test_validate_analysis_params_empty_mismatch_disables_not_rejected():
    """An explicit empty mismatch list — what `oligo.parse_mismatch_positions`
    returns for the CLI's `--mismatch 0` — is accepted as "no mismatch",
    distinct from an *unset* `mismatch` (which still defaults to `[1]`)."""
    params = validate_analysis_params(mismatch=[])
    assert params.stem_mismatch == []


@pytest.mark.parametrize("top_n", [0, 1, 100])
def test_validate_top_n_accepts_nonnegative_values(top_n):
    assert validate_top_n(top_n) == top_n


def test_validate_top_n_rejects_negative_values():
    with pytest.raises(AnalysisParameterError, match="top"):
        validate_top_n(-1)


def test_make_shrnas_short_sequence():
    """An mRNA shorter than the guide length has no valid window, so it
    yields an empty candidate list rather than raising."""
    mrnas = [("short", Seq("ATGC"))]
    params = AnalysisParams(length=22)
    result = make_shrnas(mrnas, params)
    assert result["short"] == []


def test_make_shrnas_multiple_mrnas():
    """Each mRNA gets its own independent candidate list, keyed by id."""
    mrnas = [
        ("a", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT")),
        ("b", Seq("GGGGGGGGGGGGGGGGGGGGGGGGG")),
    ]
    params = AnalysisParams(length=22)
    result = make_shrnas(mrnas, params)
    assert len(result["a"]) == 4
    assert len(result["b"]) == 4


def test_score_one_gc_filter(ref):
    """Guides outside the [gc_min, gc_max] window are rejected with reason
    "gc", on both the too-high and too-low side."""
    params = AnalysisParams(gc_min=30, gc_max=70)

    high_gc = "GGGGGGGGGGGGGGGGGGGGGG"
    candidate, reason = score_one(high_gc, "test", 0, 1, params, ref)
    assert candidate is None
    assert reason == "gc"

    low_gc = "ATATATATATATATATATATAT"
    candidate, reason = score_one(low_gc, "test", 0, 1, params, ref)
    assert candidate is None
    assert reason == "gc"


def test_score_one_motif_filter(ref):
    """A guide containing one of MOTIFS_TO_AVOID is rejected with reason
    "motif" when remove_motifs is enabled."""
    params = AnalysisParams(remove_motifs=True)

    has_tttt = "ATCAGCTAGCTTTTAGCTAGCT"
    candidate, reason = score_one(has_tttt, "test", 0, 1, params, ref)
    assert candidate is None
    assert reason == "motif"


def test_motif_filter_triggers_for_each_motif(ref):
    """Regression for the UGGC bug: every motif in the list must be reachable.
    Each motif is embedded in a 22-nt guide with neutral GC content."""
    params = AnalysisParams(remove_motifs=True)

    # Templates chosen to be 22 nt, ~50% GC, and only the motif as the trigger.
    # We embed each motif in the middle of a benign 22-nt scaffold.
    scaffold_pre = "ATCAGCTAGC"  # 10 nt, neutral
    for motif in MOTIFS_TO_AVOID:
        guide = (
            scaffold_pre + motif + "AGCT" * ((22 - len(scaffold_pre) - len(motif)) // 4)
        )
        # Pad to exactly 22
        guide = (guide + "ATGC")[:22]
        assert motif in guide, f"test setup wrong for {motif}"
        candidate, reason = score_one(guide, "test", 0, 1, params, ref)
        assert reason == "motif", f"motif {motif} was not filtered (got {reason})"
        assert candidate is None


def test_score_one_passes(ref):
    """A guide that clears every filter returns a Candidate (not None) whose
    antisense and start_pos echo the inputs."""
    params = AnalysisParams()

    guide = "ATCAGCTAGCTAGCTAGCTACG"
    candidate, reason = score_one(guide, "test", 0, 1, params, ref)
    assert candidate is not None
    assert reason is None
    assert candidate.antisense == guide
    assert candidate.start_pos == 0


def test_score_range(ref):
    """Every scored candidate's efficacy score falls within the range the
    Matveeva coefficient table can produce for a 22-mer."""
    params = AnalysisParams()

    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    candidates = make_shrnas(mrnas, params)

    for start_pos, guide in candidates["test"]:
        candidate, _ = score_one(guide, "test", start_pos, 1, params, ref)
        if candidate:
            assert -1.306 <= candidate.efficacy_score <= 1.336, (
                f"Efficacy score {candidate.efficacy_score} out of range"
            )


def test_process_shrnas_with_filters():
    """End-to-end run through process_shrnas with GC and motif filters on:
    the windows_scanned/unique_guides/passed stats stay internally consistent
    and results never exceed the number of windows generated."""
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    params = AnalysisParams(length=22, gc_min=30, gc_max=70, remove_motifs=True)
    candidates = make_shrnas(mrnas, params)
    results, stats = process_shrnas(candidates, mrnas, params)

    assert stats["windows_scanned"] == 4
    assert stats["unique_guides"] == 4
    assert stats["windows_scanned"] >= stats["unique_guides"]
    assert stats["passed"] == len(results)
    assert len(results) <= 4


def test_sort_by_efficacy_score():
    """Default params sort results descending by Candidate.efficacy_score."""
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    params = AnalysisParams()
    candidates = make_shrnas(mrnas, params)
    results, _ = process_shrnas(candidates, mrnas, params)

    for i in range(len(results) - 1):
        assert results[i].efficacy_score >= results[i + 1].efficacy_score


def test_sort_by_dg_score():
    """sort_by_dg=True switches the ranking key to Candidate.delta_g_score
    instead of the efficacy score."""
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    params = AnalysisParams(sort_by_dg=True)
    candidates = make_shrnas(mrnas, params)
    results, _ = process_shrnas(candidates, mrnas, params)

    for i in range(len(results) - 1):
        assert results[i].delta_g_score >= results[i + 1].delta_g_score


def test_reference_csv_loads_correct_shape():
    """coeffs has 22 positions × {A,C,G,T,dG}; dG_bins is 21 × 4; deltaG has 16 keys."""
    coeffs = load_coeffs()
    assert coeffs.shape[0] == 22
    for col in ("A", "C", "G", "T"):
        assert col in coeffs.columns
    bins = load_dG_bins()
    assert bins.shape == (21, 4)
    deltaG = load_deltaG()
    assert len(deltaG) == 16


def test_deltaG_has_all_16_dinucleotides():
    """The loader must symmetry-expand the on-disk table to all 16 RNA NN steps,
    so the engine never falls back to a derived key."""
    d = load_deltaG()
    expected = {a + b for a in "ACGU" for b in "ACGU"}
    assert set(d.keys()) == expected
    assert len(d) == 16


def test_load_reference_data_caches():
    """load_reference_data must memoize so repeat calls are cheap."""
    a = load_reference_data()
    b = load_reference_data()
    assert a is b


def test_score_one_rejects_wrong_length(ref):
    """A guide shorter than MIN_GUIDE_LENGTH raises ValueError rather than
    being scored against a coefficient table it doesn't fit."""
    params = AnalysisParams()
    with pytest.raises(ValueError, match="must be 20–22 nt"):
        score_one("ATCG", "test", 0, 1, params, ref)


def test_score_one_discards_guides_with_wildcards(ref):
    """A wildcard guide is a rejection, not an error: the sliding window is
    expected to straddle the N runs a reference sequence carries."""
    params = AnalysisParams()
    guide = "ATCAGCTAGCTAGCTAGCTANN"  # 22 nt, two degenerate N's
    assert len(guide) == 22
    candidate, reason = score_one(guide, "test", 0, 1, params, ref)
    assert candidate is None
    assert reason == "wildcard"


def test_score_one_still_raises_on_non_nucleotide_characters(ref):
    """Nothing upstream can produce these, so they mean a caller passed
    something wrong — unlike an IUPAC code, which is legitimate input."""
    params = AnalysisParams()
    guide = "ATCAGCTAGCTAGCTAGCTA?!"  # 22 nt, two non-nucleotides
    assert len(guide) == 22
    with pytest.raises(ValueError, match="unsupported base"):
        score_one(guide, "test", 0, 1, params, ref)


def test_score_one_gc_boundary_inclusive(ref):
    """Both GC comparisons in `score_one` are strict (`<`/`>`), so a guide
    whose GC% exactly equals gc_min or gc_max must still pass."""
    guide = "GGGGGGGGGAAAAAAAAAAAAA"  # 22 nt, 9 G/C
    assert len(guide) == 22
    gc_pct = gc_fraction(guide) * 100.0
    params = AnalysisParams(gc_min=gc_pct, gc_max=gc_pct)
    candidate, reason = score_one(guide, "test", 0, 1, params, ref)
    assert candidate is not None
    assert reason is None


def test_process_shrnas_efficacy_score_filtered_stat():
    """An unreachable min_efficacy_score rejects every candidate, and all of them are
    attributed to the "efficacy_score_filtered" stat rather than passing silently."""
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT"))]
    params = AnalysisParams(min_efficacy_score=999.0)  # unreachable: everything rejected
    candidates = make_shrnas(mrnas, params)
    results, stats = process_shrnas(candidates, mrnas, params)

    assert results == []
    assert stats["efficacy_score_filtered"] == stats["unique_guides"]
    assert stats["passed"] == 0


def test_process_shrnas_target_filtered_stat():
    """A guide unique to one mRNA must be rejected, and counted, once
    `multitarget` requires it to occur in more mRNAs than it does."""
    mrnas = [
        ("a", Seq("ATGCGTAGCTAGCTAGCTAGCTGAT")),
        ("b", Seq("GGGGGGGGGGGGGGGGGGGGGGGGG")),
    ]
    params = AnalysisParams(multitarget=2)
    candidates = make_shrnas(mrnas, params)
    results, stats = process_shrnas(candidates, mrnas, params)

    assert all(c.mrna_id != "a" for c in results)
    assert stats["target_filtered"] > 0


def test_process_shrnas_keeps_defined_stretches_around_a_wildcard():
    """The point of tolerating wildcards: an N run costs only the windows that
    overlap it, not the whole sequence."""
    left = "ATGCGTAGCTAGCTAGCTAGCTGAT"  # 25 nt, long enough for one 22-mer
    right = "GGCATTACGTACGATCAGCATCGAT"  # 25 nt
    mrnas = [("test", Seq(f"{left}NNNNN{right}"))]
    params = AnalysisParams()
    results, stats = process_shrnas(make_shrnas(mrnas, params), mrnas, params)

    assert results, "defined flanks either side of the N run must still yield guides"
    assert stats["wildcard_filtered"] > 0
    assert stats["passed"] + stats["wildcard_filtered"] <= stats["unique_guides"]
    assert all("N" not in c.antisense for c in results)


def test_process_shrnas_finds_nothing_when_every_defined_stretch_is_too_short():
    """Defined runs shorter than the guide length leave no window to read."""
    mrnas = [("test", Seq("ATGCGTAGCNATGCGTAGCNATGCGTAGCN"))]
    params = AnalysisParams()
    results, stats = process_shrnas(make_shrnas(mrnas, params), mrnas, params)

    assert results == []
    assert stats["wildcard_filtered"] == stats["unique_guides"]
    assert stats["passed"] == 0


def test_score_one_no_keyerror_on_arbitrary_guides(ref):
    """Every 22-mer of A/C/G/T must score without KeyError after the
    deltaG symmetry expansion (regression for the dinuc-fallback bug)."""
    params = AnalysisParams()

    guides = [
        "AAAAAAAAAAAAAAAAAAAAAA",
        "CCCCCCCCCCCCCCCCCCCCCC",
        "GGGGGGGGGGGGGGGGGGGGGG",
        "TTTTTTTTTTTTTTTTTTTTTT",
        "ACACACACACACACACACACAC",
        "AGAGAGAGAGAGAGAGAGAGAG",
        "TGTGTGTGTGTGTGTGTGTGTG",
        "ATGCATGCATGCATGCATGCAT",
    ]
    for g in guides:
        # Should not raise KeyError
        score_one(g, "test", 0, 1, params, ref)


def test_make_shrnas_rejects_duplicate_mrna_identifiers():
    """Candidates are keyed by identifier, so a duplicate would overwrite an
    earlier entry and silently drop a whole sequence's candidates. The input
    layer rejects duplicates first; this is the library-API backstop."""
    mrnas = [
        ("gene", Seq("ATGCATGCATGCATGCATGCATGC")),
        ("gene", Seq("GGCCGGCCGGCCGGCCGGCCGGCC")),
    ]
    with pytest.raises(ValueError, match="Duplicate mRNA identifier"):
        make_shrnas(mrnas, AnalysisParams())


@pytest.mark.parametrize("length", [20, 21, 22])
def test_score_one_accepts_every_supported_guide_length(ref, length):
    """Every length in MIN..MAX scores without error. The coefficient table is
    5'-anchored, so a short guide uses the first `length` rows."""
    guide = "ATGCATGCATGCATGCATGCAT"[:length]
    candidate, reason = score_one(
        guide=guide,
        mrna_id="m",
        start_pos=0,
        candidate_idx=1,
        params=AnalysisParams(length=length),
        ref=ref,
    )
    assert reason is None
    assert candidate is not None
    assert len(candidate.antisense) == length


@pytest.mark.parametrize("length", [19, 23])
def test_score_one_rejects_lengths_outside_the_supported_range(ref, length):
    """Outside MIN..MAX there is no coefficient row to score against (too long)
    or the model is being stretched past what it covers (too short)."""
    guide = ("ATGCATGCATGCATGCATGCATGC" * 2)[:length]
    with pytest.raises(ValueError, match="must be 20–22 nt"):
        score_one(
            guide=guide,
            mrna_id="m",
            start_pos=0,
            candidate_idx=1,
            params=AnalysisParams(length=length),
            ref=ref,
        )


def test_shorter_guides_score_lower_because_the_score_is_a_positional_sum():
    """Guard for the documented caveat: the efficacy score sums per-position
    terms, so it is NOT comparable across guide lengths. A truncated guide has
    fewer terms. Anyone who "fixes" this into a mean must update the docs."""
    ref_data = load_reference_data()
    full = "ATGCATGCATGCATGCATGCAT"
    short = full[:20]
    scored_full, _ = score_one(
        guide=full,
        mrna_id="m",
        start_pos=0,
        candidate_idx=1,
        params=AnalysisParams(length=22),
        ref=ref_data,
    )
    scored_short, _ = score_one(
        guide=short,
        mrna_id="m",
        start_pos=0,
        candidate_idx=1,
        params=AnalysisParams(length=20),
        ref=ref_data,
    )
    dropped = sum(ref_data.coeffs.loc[pos, full[pos]] for pos in range(20, 22))
    assert scored_full.efficacy_score == pytest.approx(
        scored_short.efficacy_score + dropped, abs=1e-4
    )
