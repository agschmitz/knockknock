# Architecture — how the algorithm actually works

Read this before touching `engine.py`, `oligo.py`, `cutsite.py`,
`backbone.py`, or `reference.py`. For *usage* (CLI flags, TUI key
bindings), see [`CLI_MANUAL.md`](CLI_MANUAL.md) and
[`TUI_MANUAL.md`](TUI_MANUAL.md) instead — this doc is about the math and
the data model, not the interface.

## 1. Candidate generation and scoring (`engine.py`)

For each input mRNA (`Bio.Seq` from a parsed FASTA record):

1. **`make_shrnas`** takes the reverse complement of the mRNA and slides a
   window of `AnalysisParams.length` (default 22 nt) across it, one base at
   a time. Each window is a candidate guide (antisense) sequence, paired
   with its 0-based start position on the reverse complement.
2. **`process_shrnas`** dedupes identical guide sequences across all input
   mRNAs (`seen_guides`), then calls **`score_one`** on each unique guide:
   - **Wildcard filter** — reject if the window contains a degenerate IUPAC
     code (`N`, `R`, `Y`, ...). `seqinput.load_mrnas` lets wildcards
     through the input layer on purpose — rejecting the record would throw
     away its well-defined regions too — so the sliding window is expected
     to produce some windows straddling an `N` run, and this is where they
     are dropped. Only the guides that actually touch an undefined base are
     lost; a sequence whose every defined stretch is shorter than the guide
     length yields nothing, which is the intended outcome. A character that
     is not a nucleotide at all still raises, since nothing upstream can
     produce one.
   - **GC filter** — reject if `gc_fraction(guide) * 100` falls outside
     `[gc_min, gc_max]`.
   - **Efficacy score** — sum `ref.coeffs.loc[position, base]` over every
     position of the guide (Matveeva et al. 2012 position-specific
     efficacy coefficients). Reject if below `min_efficacy_score` (when set
     — `None` means no floor, which is the `lax`/`standard` preset default).
   - **Motif filter** (optional, `remove_motifs`) — reject if any of
     `MOTIFS_TO_AVOID = ("TTTT", "AAAA", "GGGG", "TGGC", "CGGA")` appears in
     the guide (DNA alphabet — the guide is checked before RNA
     transcription).
   - **Duplex ΔG** — transcribe the guide to RNA, walk overlapping
     dinucleotides, sum `ref.deltaG[dinuc]` (Xia et al. 1998
     nearest-neighbor free energies) plus the `DG_INITIATION = 4.09`
     kcal/mol correction → `dg`. In parallel, bin each dinucleotide's ΔG
     into one of 4 stability classes via `_score_bins` and sum the
     matching `ref.dG_bins` cell → `delta_g_score` (a *separate*, ΔG score —
     not the same thing as the raw `dg` value).
   - **Multi-target filter** — applied after scoring, across all mRNAs: a
     guide's `targets` count is how many of the input mRNAs' reverse
     complements contain it as a substring (checked against a
     precomputed set per mRNA, `_build_revcomp_lookup`, so this is O(1)
     per mRNA instead of a fresh substring search each time). Reject if
     `targets < params.multitarget`.
   - The passenger strand for each surviving candidate is computed
     immediately via `oligo.compute_mismatched_passenger` (see §2) and
     stored as `Candidate.sense`, even before the user reaches the oligo
     designer — so exports always have a sense-strand column.
3. Results are sorted descending by `efficacy_score` (default) or
   `delta_g_score` (`sort_by_dg=True`), and `process_shrnas` returns
   `(candidates, stats)` where `stats` counts `windows_scanned`,
   `unique_guides`, and rejections by reason (`wildcard_filtered`,
   `gc_filtered`, `efficacy_score_filtered`, `motif_filtered`,
   `target_filtered`), plus `passed`. Those counts are what
   `export.no_candidates_message` reads to explain an empty result — no
   output path prints a bare header row.

`QUALITY_PRESETS` (`lax` / `standard` / `strict`) are just canned
`(gc_min, gc_max, remove_motifs, min_efficacy_score)` bundles applied in
`cli.py` before individual `--gc-min`-style flags override them.

Reference tables are loaded once per process via
`reference.load_reference_data()` (`@lru_cache(maxsize=1)`) into a
`ReferenceData(coeffs, deltaG, dG_bins)` bundle. Column-level detail
(shapes, value ranges, source citations) lives in
[`knockknock/data/README.md`](../knockknock/data/README.md) — not
duplicated here.

## 2. Oligo assembly (`oligo.py`, `cutsite.py`, `backbone.py`)

This subsystem turns one guide + a backbone into two synthesizable oligo
strings. It's the part of the codebase with the most non-obvious
invariants, because restriction-site geometry is easy to get backwards.

### 2.1 Cut sites and overhangs (`cutsite.py`)

A backbone junction is specified in conventional notation — the
recognition sequence with a slash at the **top-strand** cut position,
e.g. `G/AATTC` (EcoRI) or `CTGCA/G` (PstI). Because Type II enzymes cut
symmetrically about the center of their recognition site, that single
slash fully determines the bottom-strand cut too: with site length `L`
and `c` bases before the slash, `k = L - 2c` gives a 5' overhang of `k`
nt (`k > 0`), blunt (`k == 0`), or a 3' overhang of `|k|` nt (`k < 0`).

Type IIS enzymes (BsaI, BbsI, BsmBI, SapI) cut *outside* their
recognition site and can't be expressed this way — their overhang comes
from the plasmid's spacer, so it's stated directly instead: `5'AATT`,
`3'TGCA`, or `blunt`. `parse_overhang` accepts either form.

**Degenerate codes are rejected only where they land on the overhang.**
`N/TCGAG` and ApoI's `R/AATTY` still pin the overhang down exactly
(`TCGA`, `AATT`), so they parse; DraIII's `CACNNN/GTG` does not, because
there the N region *is* the overhang and no oligo can be ordered against
it. `validate_dna_sequence(..., allow_degenerate=True)` covers the
recognition sequence, `Overhang.__post_init__` enforces the stricter rule
on the overhang alone. A partly-degenerate site is still searched for
correctly inside the insert — `site_pattern`/`site_occurs_in` expand
IUPAC codes to character classes rather than matching literally.

**A site always describes the digested backbone, never the insert.** A
directly-stated overhang is the backbone's own single-stranded end,
written 5'→3' along whichever strand protrudes there — the way an
overhang is written anywhere else. The insert is then designed to fit
it, so the typed sequence turns up on the oligo only as its reverse
complement. Which strand of the backbone protrudes flips between the
two junctions:

| junction | kind | backbone strand that protrudes | stated 5'→3' runs |
|---|---|---|---|
| left (`--site-5`) | 5' | bottom | right-to-left |
| left (`--site-5`) | 3' | top | left-to-right |
| right (`--site-3`) | 5' | top | left-to-right |
| right (`--site-3`) | 3' | bottom | right-to-left |

That is why `parse_overhang` takes a `junction` (`LEFT`/`RIGHT`) and why
`Overhang` carries one: the same text means two different sequences at
the two ends. Cut-site notation is exempt — `G/AATTC` names the *uncut
plasmid's* top strand, which is already top-strand coordinates wherever
the junction sits.

Every real Type II overhang is palindromic (the recognition site is), so
this only ever bites on stated Type IIS / synthetic ends — which is
exactly how it went unnoticed until a non-palindromic `5'TTAT` was
entered.

**Everything is stored in top-strand coordinates.** An `Overhang.top` is
always the sequence read 5'→3' on the top strand, even when the
physically-protruding strand is the bottom one. This is what lets
`build_strands`/`make_duplex` treat "draw the top strand left-to-right"
as "emit `.top` verbatim" uniformly, regardless of which end of the
insert or which physical strand is involved — don't reintroduce a
"raw overhang string" representation that loses this distinction.
`Overhang.backbone_top` is the counterpart for display: the same end
restated as the backbone's 5'→3' overhang, which is what `describe()`
shows and what a site field round-trips.

### 2.2 Assembling the insert (`oligo.py`)

`OligoConfig.core(arm1, arm2)` builds the fully double-stranded body:
`5'_flank + arm1 + loop + arm2 + 3'_flank`, where `arm1`/`arm2` are
guide/passenger in the order `guide_on_5p` selects.

**`keep_sites` brackets that body with the rest of each recognition
sequence.** Digestion splits the site across the two fragments: the
overhang is shared, but the bases outside it are not. The left fragment
keeps `Overhang.site_head`, the right keeps `Overhang.site_tail`, and the
insert *is* the right fragment at its own 5' junction and the left one at
its 3' junction — so `Overhang.insert_site_bases` returns `site_tail` at
a LEFT junction and `site_head` at a RIGHT one. Carrying them back is
what makes ligation restore the site instead of scarring it: XhoI's
`C/TCGAG` at a 5' junction gives a forward oligo of `TCGA` + `G` + core,
and the backbone's own `C` completes `CTCGAG`.

These bases are double-stranded, so they belong to the core rather than
the overhang — which is why they extend the annealed overlap, the Tm, and
both oligos, with only the overhangs left single-stranded. The setting is
stored per backbone (`backbones.csv`'s `keep_sites` column, on only for
pINDUCER), since whether a vector's sites are normally regenerated is a
property of that vector; `OligoConfig.from_backbone(keep_sites=None)`
takes the stored value and `--keep-sites`/`--no-keep-sites` override it
per run. That tri-state matters: a plain `False` default would silently
switch the flag off for every backbone that stores `True`.
`check_sites` reports the two ends where the flag cannot deliver: a directly stated overhang has no site to rebuild, and a
degenerate base in the insert's own contribution (`N/TCGAG` at a *3'*
junction, where the insert would owe the `N`) makes it unorderable.

`build_strands` then prepends/appends the correct overhang to the
correct oligo, per this table (also in the `build_strands` docstring):
a 5'-protruding left junction extends the **forward** oligo's 5' end; a
3'-protruding left junction extends the **reverse** oligo's 3' end (via
reverse-complementing the top-strand overhang); the right junction is
the mirror image. Blunt junctions contribute nothing.

`make_duplex` lays the same insert out on a shared column grid
(`InsertDuplex`) — this is a pure rendering helper, kept separate from
`build_strands` so the two never drift on which junction contributes to
which oligo. `duplex_view.py` builds on `InsertDuplex` to draw the TUI's
color-coded four-row junction diagram (`render_diagram`); it only takes
`OligoConfig` and the two arms, so it stays a pure function too — no
Textual import, testable and reusable outside the TUI.

### 2.3 The passenger-strand mismatch

`compute_mismatched_passenger(guide, mismatch_positions, guide_on_5p=False)`
reverse-complements the guide, then at each 1-based `mismatch_positions`
entry replaces the base via the deterministic transversion table
`MISMATCH_TABLE = {"G": "T", "T": "G", "A": "C", "C": "A", "U": "G"}`.
This is intentionally **not random** — an earlier implementation chose
the substitute base randomly, which made re-running the tool on the same
input produce different oligos. Keep it deterministic.

Positions are counted **from the base of the stem** (the end nearest the
flanking/backbone sequence, farthest from the loop) — position 1 always
opens up the base of the stem, never the loop end. Which end of the
*passenger string itself* that is depends on `guide_on_5p`: `OligoConfig
.core` always builds `arm1` with its front adjacent to the 5' flank (the
base of the stem) and `arm2` with its front adjacent to the loop (its
*back* is the base of the stem, next to the 3' flank). So when the guide
sits on arm2 (`guide_on_5p=False`, the default — the passenger is arm1),
position 1 indexes the passenger from its own front; when the guide sits
on arm1 (`guide_on_5p=True` — the passenger is arm2), position 1 must
index the passenger from its own *back* instead, or it would land at the
loop end. `assemble_oligos`/`assemble_batch` pass `config.guide_on_5p`
through for this reason. `engine.score_one` calls this function at the
default `guide_on_5p=False` (no backbone/arm order is chosen yet at
scoring time, so there is no stem to speak of) purely to fill
`Candidate.sense` with a preview passenger.

`duplex_view.core_parts` mirrors the same base-of-stem convention for the
diagram's mismatch highlight, independently of which content occupies
each arm: `mismatch_positions` apply as-is (front-counted) to arm1 and
mirrored (`len(arm) - p + 1`) to arm2, so the highlighted column is always
the actual base of the stem on both drawn arms, for either orientation.

### 2.4 Tm

`_annealed_tm` computes the melting temperature via
`Bio.SeqUtils.MeltingTemp.Tm_NN` (SantaLucia 1998 nearest-neighbor
defaults) over the **core only** — i.e. excluding both overhangs, which
stay single-stranded regardless of which end protrudes. It deliberately
ignores the seed-region mismatch's effect on Tm (small and consistent
across candidates, per the `oligo.py` docstring).

### 2.5 Validation warnings

`validate_design` = `check_ends` (are the two junctions interchangeable,
i.e. is cloning non-directional — true whenever both ends are blunt, or
both are the same kind with matching sequence) + `check_sites` (where
`keep_sites` could not rebuild the recognition sequence) + `check_insert`
(does either backbone enzyme's recognition site also turn up inside the
assembled insert). All return warning strings; an empty list means no
issues. These are advisory — nothing here blocks assembly.

`check_sites` reports only *failures*: a junction rebuilt cleanly says
nothing, since the restored bases are visible in the oligo and coloured
in the diagram. A junction that silently could not be rebuilt is the
dangerous case, so that one always speaks up.

### 2.6 Backbones (`backbone.py`)

`load_backbones()` reads `knockknock/data/backbones.csv` once (module-
level cache, reset by `save_backbone`), parsing each row's junction
columns through `cutsite.parse_overhang` — `5_prime_site` as `LEFT`,
`3_prime_site` as `RIGHT`, since both columns describe the backbone's
own ends (see 2.1). Boolean columns (`guide_on_5p`, `keep_sites`) go
through `_truthy(row.get(...))`: `.get`, not `[...]`, because the user's
own database under `~/.config` is seeded once and then only appended to,
so an installed copy can predate any column added later and must keep
loading rather than `KeyError`. `save_backbone` appends a new row and invalidates
the cache; it refuses to overwrite an existing name. The row stores
`Overhang.spec`, the text as typed, so a saved backbone re-parses to the
same junction it was created with.

## 3. Export (`export.py`, `oligo.py`'s `export_oligos_*`)

`export.py` covers the scoring-command outputs (terminal table, CSV, XLSX,
JSON, FASTA) plus provenance metadata (`_metadata`: version, timestamp,
input SHA-256, and every `AnalysisParams` field) written into a CSV
comment header, a JSON `meta` key, or an XLSX "Run Info" sheet. XLSX
conditional formatting uses `xlsxwriter`-specific APIs (column width,
3-color scales) — swapping the pandas Excel engine would break it
silently.

`to_terminal`/`to_fasta`/`to_json_text`/`to_xlsx` each take an optional
`oligo_pairs` (and `backbone_name`), used when the CLI's `--backbone`/manual
insert flags auto-trigger oligo design alongside scoring. Each format keeps
the oligo data as its own section rather than merging it into the candidate
rows: a separate panel in the terminal table, extra FASTA records, an
`"oligos"` key alongside `"candidates"` in JSON, or a second XLSX sheet.
CSV is the exception — it can't hold two tables in one file, so the CLI
writes a second `*_oligos.csv` next to the candidates CSV instead of
passing anything through `to_csv_text`.

Oligo-pair export (`oligo.oligo_rows` and its CSV/JSON/text wrappers) is
kept in `oligo.py`, not `export.py`, and both the CLI and the TUI call
the same `oligo_rows` so the two paths can't drift on which fields an
oligo export carries.
