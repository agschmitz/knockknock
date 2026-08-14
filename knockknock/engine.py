"""Core pipeline: generate → score → filter → sort."""

from dataclasses import dataclass, field
from math import isfinite

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from knockknock.cutsite import invalid_dna_chars
from knockknock.oligo import compute_mismatched_passenger
from knockknock.reference import ReferenceData, load_reference_data

MOTIFS_TO_AVOID: tuple[str, ...] = ("TTTT", "AAAA", "GGGG", "TGGC", "CGGA")
DG_INITIATION = 4.09

#: Guide lengths the bundled reference tables can score.
#:
#: The Matveeva coefficient table has 22 rows and the ΔG bin table 21, both
#: indexed from the guide's 5' end. A shorter guide therefore uses the first
#: ``n`` rows — the model is **5'-anchored**: position 1 of a 20-mer is scored
#: with the same coefficients as position 1 of a 22-mer, and the 3'-most rows
#: simply go unused.
#:
#: Consequence worth knowing before comparing results: the efficacy score is a
#: *sum* over positions, so it is not comparable across lengths. A 20-mer has
#: two fewer terms and will tend to score lower than an equally good 22-mer.
#: Rank candidates within one length, not across a mixed-length run.
# The design workflow deliberately excludes 19-mers.  The bundled scoring
# tables can technically score them, but discovery is restricted to 20–22 nt
# guides so every generated target has the intended minimum length.
MIN_GUIDE_LENGTH = 20
MAX_GUIDE_LENGTH = 22
DEFAULT_GUIDE_LENGTH = 22

#: Backwards-compatible alias for the previous single supported length.
SUPPORTED_GUIDE_LENGTH = DEFAULT_GUIDE_LENGTH


QUALITY_PRESETS: dict[str, dict] = {
    "lax": {
        "gc_min": 0.0,
        "gc_max": 100.0,
        "remove_motifs": False,
        "min_efficacy_score": None,
    },
    "standard": {
        "gc_min": 20.0,
        "gc_max": 80.0,
        "remove_motifs": True,
        "min_efficacy_score": None,
    },
    "strict": {
        "gc_min": 30.0,
        "gc_max": 60.0,
        "remove_motifs": True,
        "min_efficacy_score": 0.0,
    },
}


@dataclass
class Candidate:
    name: str
    #: 0-based offset into the **reverse complement** of the source mRNA
    #: (i.e. `make_shrnas`' sliding-window coordinate), not the original
    #: input sequence.
    start_pos: int
    antisense: str
    sense: str
    efficacy_score: float
    delta_g_score: float
    gc_pct: float
    dg: float
    targets: int
    mrna_id: str


@dataclass
class AnalysisParams:
    length: int = DEFAULT_GUIDE_LENGTH
    gc_min: float = 0.0
    gc_max: float = 100.0
    min_efficacy_score: float | None = None
    sort_by_dg: bool = False
    remove_motifs: bool = False
    multitarget: int = 1
    rna_output: bool = False
    stem_mismatch: list[int] = field(default_factory=lambda: [1])
    # Reserved for future stochastic features. Currently no engine code
    # consumes randomness, but the value is recorded in metadata exports
    # (Phase 5.2) so reproducible reruns are possible later.
    seed: int | None = None


class AnalysisParameterError(ValueError):
    """An analysis setting is outside the range accepted by the engine."""


def validate_analysis_params(
    *,
    length: int = DEFAULT_GUIDE_LENGTH,
    gc_min: float = 0.0,
    gc_max: float = 100.0,
    min_efficacy_score: float | None = None,
    sort_by_dg: bool = False,
    remove_motifs: bool = False,
    multitarget: int = 1,
    rna_output: bool = False,
    mismatch: list[int] | None = None,
    seed: int | None = None,
) -> AnalysisParams:
    """Validate shared analysis settings and return engine-ready parameters.

    The CLI and TUI both use this function.  Keeping the rules at the engine
    boundary prevents the two interfaces from accepting different filters.

    ``mismatch`` is a list of 1-based passenger-strand mismatch positions
    (default ``[1]``) rather than a single position, so the same value the
    CLI parses for scoring can be reused unchanged for oligo design. An
    explicit empty list (from `oligo.parse_mismatch_positions` parsing
    ``"0"``) disables mismatch introduction entirely rather than erroring —
    only an *unset* ``mismatch`` falls back to the ``[1]`` default.
    """
    if not MIN_GUIDE_LENGTH <= length <= MAX_GUIDE_LENGTH:
        raise AnalysisParameterError(
            f"length must be between {MIN_GUIDE_LENGTH} and {MAX_GUIDE_LENGTH}; "
            "the bundled scoring coefficients only cover that range."
        )
    if not isfinite(gc_min) or not isfinite(gc_max):
        raise AnalysisParameterError("GC limits must be finite numbers.")
    if not (0 <= gc_min < gc_max <= 100):
        raise AnalysisParameterError(
            "GC limits must satisfy 0 ≤ minimum < maximum ≤ 100."
        )
    if min_efficacy_score is not None and not isfinite(min_efficacy_score):
        raise AnalysisParameterError("Minimum efficacy score must be a finite number.")
    if multitarget < 1:
        raise AnalysisParameterError("multitarget must be at least 1.")
    stem_mismatch = list(mismatch) if mismatch is not None else [1]
    bad = [m for m in stem_mismatch if not 1 <= m <= length]
    if bad:
        raise AnalysisParameterError(
            f"mismatch position(s) {bad} must be between 1 and {length} "
            "(pass an empty list to disable mismatch introduction)."
        )
    return AnalysisParams(
        length=length,
        gc_min=gc_min,
        gc_max=gc_max,
        min_efficacy_score=min_efficacy_score,
        sort_by_dg=sort_by_dg,
        remove_motifs=remove_motifs,
        multitarget=multitarget,
        rna_output=rna_output,
        stem_mismatch=stem_mismatch,
        seed=seed,
    )


def validate_top_n(top_n: int) -> int:
    """Validate the shared export limit."""
    if top_n < 0:
        raise AnalysisParameterError("top must be zero or greater.")
    return top_n


def _score_bins(position: int, dG_NN: float, dG_bins: pd.DataFrame) -> float:
    if position >= dG_bins.shape[0]:
        return 0.0
    if dG_NN > -1:
        return dG_bins.iloc[position, 0]
    elif dG_NN > -2:
        return dG_bins.iloc[position, 1]
    elif dG_NN > -3:
        return dG_bins.iloc[position, 2]
    else:
        return dG_bins.iloc[position, 3]


def make_shrnas(
    mrnas: list[tuple[str, Seq]], params: AnalysisParams
) -> dict[str, list[tuple[int, str]]]:
    """Slide a `params.length`-nt window across each mRNA's reverse complement.

    Returns, per mRNA id, a list of ``(start_pos, guide)`` pairs — every
    possible guide (antisense) candidate. ``start_pos`` is 0-based **on the
    reverse complement**, not on the original input sequence.

    Raises ``ValueError`` on a repeated mRNA id. The return value is keyed by
    id, so a duplicate would overwrite the earlier entry and silently discard
    a whole input sequence's candidates. `seqinput.load_mrnas` already rejects
    duplicates at the input boundary; this is the backstop for direct callers
    of the library API.
    """
    result: dict[str, list[tuple[int, str]]] = {}
    for mrna_id, seq in mrnas:
        if mrna_id in result:
            raise ValueError(
                f"Duplicate mRNA identifier {mrna_id!r}: candidates are keyed "
                "by identifier, so every input sequence needs a unique one."
            )
        rc = seq.upper().reverse_complement()
        n_windows = len(rc) - params.length + 1
        if n_windows < 1:
            result[mrna_id] = []
            continue
        candidates = [
            (pos, str(rc[pos : pos + params.length])) for pos in range(n_windows)
        ]
        result[mrna_id] = candidates
    return result


def _build_revcomp_lookup(mrnas: list[tuple[str, Seq]], length: int) -> list[set[str]]:
    sets: list[set[str]] = []
    for _, seq in mrnas:
        rc = str(seq.upper().reverse_complement())
        n = len(rc) - length + 1
        if n < 1:
            sets.append(set())
        else:
            sets.append({rc[i : i + length] for i in range(n)})
    return sets


def candidate_sort_key(candidate: Candidate, sort_by_dg: bool) -> float:
    """Rank value for one candidate: `delta_g_score` when sorting by duplex
    stability, else the `efficacy_score`.

    Shared by `process_shrnas` and the TUI's live filter bar (`tui.
    ExploreApp._filtered`) so a `knockknock run --sort-dg` and pressing `s`
    in the explorer order candidates identically.
    """
    return candidate.delta_g_score if sort_by_dg else candidate.efficacy_score


def rejection_reason(
    gc_pct: float,
    efficacy_score: float,
    guide: str,
    params: AnalysisParams,
    *,
    targets: int | None = None,
) -> str | None:
    """The name of the first filter that rejects these values, else ``None``.

    **The single definition of what the GC / efficacy score / motif /
    multi-target filters mean.** `score_one` calls it while scoring, and the
    TUI's live filter bar calls it against already-scored candidates — so the
    interactive view and a `knockknock run` with the same settings can never
    disagree about which candidates pass.

    ``targets`` is optional because `score_one` does not know it yet: the
    multi-target count needs every mRNA's candidates at once, so
    `process_shrnas` supplies it on a second call. Pass ``None`` to skip that
    filter.
    """
    if gc_pct < params.gc_min or gc_pct > params.gc_max:
        return "gc"
    if params.min_efficacy_score is not None and efficacy_score < params.min_efficacy_score:
        return "efficacy_score"
    if params.remove_motifs and any(m in guide for m in MOTIFS_TO_AVOID):
        return "motif"
    if targets is not None and targets < params.multitarget:
        return "target"
    return None


def score_one(
    guide: str,
    mrna_id: str,
    start_pos: int,
    candidate_idx: int,
    params: AnalysisParams,
    ref: ReferenceData,
) -> tuple[Candidate | None, str | None]:
    """Score one guide and apply the wildcard/GC/efficacy-score/motif filters,
    in order.

    Returns ``(candidate, None)`` on success, or ``(None, reason)`` where
    ``reason`` is one of ``"wildcard"``, ``"gc"``, ``"efficacy_score"``, or
    ``"motif"`` — whichever filter rejected it first. The multi-target filter
    is not applied here; it needs every mRNA's candidates at once, so
    `process_shrnas` applies it afterward.

    A guide containing a degenerate IUPAC code is a **rejection, not an
    error**: `seqinput.load_mrnas` deliberately lets wildcards through so the
    defined stretches around them stay targetable, which means the sliding
    window is expected to produce some windows that straddle one. Characters
    that are not nucleotides at all still raise — nothing upstream can produce
    those, so they mean a caller passed something wrong.
    """
    if not MIN_GUIDE_LENGTH <= len(guide) <= MAX_GUIDE_LENGTH:
        raise ValueError(
            f"Guide length must be {MIN_GUIDE_LENGTH}–{MAX_GUIDE_LENGTH} nt; "
            f"got {len(guide)}."
        )
    degenerate, other = invalid_dna_chars(guide.upper())
    if other:
        raise ValueError(
            f"Guide contains unsupported base(s) {other!r}; "
            "expected DNA bases A, C, G, T."
        )
    if degenerate:
        return None, "wildcard"
    guide = guide.upper()
    guide_seq = Seq(guide)

    gc_pct = gc_fraction(guide_seq) * 100.0

    efficacy_score = 0.0
    for pos in range(len(guide)):
        nt = guide[pos]
        efficacy_score += ref.coeffs.loc[pos, nt]

    # GC / min-efficacy-score / motif, in that order — see `rejection_reason`,
    # which the TUI's filter bar shares so the two views can't disagree.
    reason = rejection_reason(gc_pct, efficacy_score, guide, params)
    if reason is not None:
        return None, reason

    # Dinucleotide ΔG scoring
    guide_rna = str(guide_seq.transcribe())
    dg = 0.0
    delta_g_score = 0.0
    for dinuc_pos in range(len(guide_rna) - 1):
        dinuc = guide_rna[dinuc_pos : dinuc_pos + 2]
        dG_NN = ref.deltaG[dinuc]
        dg += dG_NN
        delta_g_score += _score_bins(dinuc_pos, dG_NN, ref.dG_bins)
    dg += DG_INITIATION

    passenger = compute_mismatched_passenger(str(guide_seq), params.stem_mismatch)
    name = f"{mrna_id}-sh{candidate_idx}"

    return Candidate(
        name=name,
        start_pos=start_pos,
        antisense=str(guide_seq),
        sense=passenger,
        efficacy_score=round(efficacy_score, 4),
        delta_g_score=round(delta_g_score, 1),
        gc_pct=round(gc_pct, 1),
        dg=round(dg, 1),
        targets=0,
        mrna_id=mrna_id,
    ), None


def process_shrnas(
    candidates_by_mrna: dict[str, list[tuple[int, str]]],
    all_mrnas: list[tuple[str, Seq]],
    params: AnalysisParams,
) -> tuple[list[Candidate], dict[str, int]]:
    """Dedupe, score, and multi-target-filter every candidate, then rank.

    Returns ``(candidates, stats)``: ``candidates`` is sorted descending by
    `Candidate.efficacy_score` (or `Candidate.delta_g_score` when
    `params.sort_by_dg`). ``stats`` counts ``windows_scanned`` (every window,
    including duplicate guides), ``unique_guides`` (post-dedup), each
    rejection reason (``wildcard_filtered``, ``gc_filtered``,
    ``efficacy_score_filtered``, ``motif_filtered``, ``target_filtered``), and
    ``passed``.
    """
    ref = load_reference_data()

    revcomp_lookup = _build_revcomp_lookup(all_mrnas, params.length)

    stats: dict[str, int] = {
        "windows_scanned": 0,
        "unique_guides": 0,
        "wildcard_filtered": 0,
        "gc_filtered": 0,
        "efficacy_score_filtered": 0,
        "motif_filtered": 0,
        "target_filtered": 0,
        "passed": 0,
    }

    scored: list[Candidate] = []
    seen_guides: set[str] = set()
    for mrna_id, candidates in candidates_by_mrna.items():
        for sh_idx, (start_pos, guide) in enumerate(candidates, start=1):
            stats["windows_scanned"] += 1
            if guide in seen_guides:
                continue
            seen_guides.add(guide)
            stats["unique_guides"] += 1
            candidate, reject_reason = score_one(
                guide=guide,
                mrna_id=mrna_id,
                start_pos=start_pos,
                candidate_idx=sh_idx,
                params=params,
                ref=ref,
            )

            if candidate is None:
                if reject_reason == "wildcard":
                    stats["wildcard_filtered"] += 1
                elif reject_reason == "gc":
                    stats["gc_filtered"] += 1
                elif reject_reason == "efficacy_score":
                    stats["efficacy_score_filtered"] += 1
                elif reject_reason == "motif":
                    stats["motif_filtered"] += 1
                continue

            # Multi-target filter: deferred to here because it needs every
            # mRNA's windows at once, which `score_one` doesn't have.
            targets = sum(1 for rc_set in revcomp_lookup if guide in rc_set)
            if (
                rejection_reason(
                    candidate.gc_pct,
                    candidate.efficacy_score,
                    guide,
                    params,
                    targets=targets,
                )
                == "target"
            ):
                stats["target_filtered"] += 1
                continue

            candidate.targets = targets
            scored.append(candidate)

    stats["passed"] = len(scored)

    scored.sort(key=lambda c: candidate_sort_key(c, params.sort_by_dg), reverse=True)

    return scored, stats
