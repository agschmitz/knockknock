# KnockKnock — CLI Manual

Full reference for the `knockknock` command line. For the interactive terminal screens, see [TUI_MANUAL.md](TUI_MANUAL.md).
## The shape of a command

```bash
knockknock [INPUT]... [OPTIONS]
```

`knockknock` has a short alias, `kk`. There are no subcommands; the first token is the input.

`INPUT` is a path to a plain-text FASTA file, or to a text file holding one bare sequence. `-` reads the same kinds of text from standard input. Without `INPUT`, KnockKnock opens the [TUI](TUI_MANUAL.md) at an interactive terminal and scores stdin on a pipe.

```bash
# interactive: file picker, then explorer
knockknock

# score a file, print a table
knockknock <mRNAs.fasta>
# read a sequence or FASTA from stdin
knockknock -
# input from stdin, without the explicit '-' when used in a pipe
cat <mRNAs.fasta> | knockknock      
```

The last two are not equivalent. `knockknock -` is an ordinary scoring run with the usual defaults ([`--quality strict`](#filtering-candidates), top 10). Bare `knockknock` on a pipe applies no filters and emits every candidate as FASTA.

## Scoring mRNAs

Scoring slides a window of `-l`/`--length` nucleotides (20–22, default 22) along the reverse complement of every input mRNA and evaluates each window. Each candidate carries two scores:

- the **efficacy score**, a sum of position-specific nucleotide coefficients from Matveeva et al. (2012). Higher is better.
- the **ΔG score**, a binned sum derived from the Xia et al. (1998) nearest-neighbour duplex free energy. Higher is better. It is distinct from **ΔG**, the raw duplex energy in kcal/mol.

Both scores are sums over positions, so they are comparable only between candidates of the same guide length. Changing `-l` changes the scale.

Candidates are ranked by efficacy score; `--sort-dg` ranks by ΔG score instead. `-t`/`--top` limits how many are kept, 10 by default. See [Designing cloning oligos](#designing-cloning-oligos) for how oligo design changes that default.

When the result is empty, KnockKnock prints the reason instead of an empty table. There are three reasons: no input sequence is long enough to hold one window, every defined A/C/G/T stretch is shorter than one window, or no candidate passed the filters.

### Filtering candidates

`--quality` selects a filter preset. The default is `strict`.

| Preset | GC range | Motif filter | Minimum efficacy score |
|---|---|---|---|
| `lax` | 0–100 % | off | none |
| `standard` | 20–80 % | on | none |
| `strict` (default) | 30–60 % | on | 0.0 |

A filter flag given explicitly overrides the preset for that filter only; the rest of the preset still applies. `--quality lax --gc-min 40` therefore means the lax preset with a GC floor of 40 %. The individual flags are `--gc-min`, `--gc-max`, `-s`/`--min-efficacy-score`, and `--remove-motifs`, which rejects homopolymer runs and other low-complexity motifs.

`--multitarget N` is not part of the presets. It requires a guide to occur in at least N of the input mRNAs and defaults to all of them, so a FASTA of transcript variants or orthologues yields only the guides that hit all of them. A value above the number of input sequences is reduced to that number, with a warning.

Degenerate IUPAC codes in the input are accepted: the sequence loads with a warning, and only the windows that overlap a wildcard base are dropped. The defined stretches around them are still scored. Characters outside the IUPAC alphabet are an error.

## Output formats and files

`-f`/`--format` accepts `table`, `fasta`, `csv`, `xlsx`, and `json`. Without `-f`, the format is `table` when stdout is a terminal and `fasta` when it is piped or redirected. An explicit `-f` always applies.

| Format | Contents |
|---|---|
| `table` | Terminal table: rank, guide, both scores, GC %, ΔG, target count, position. Below it, a summary panel with the input, the active filters, and the window counts. Not written to a file. |
| `fasta` | One record per guide; rank, efficacy score, GC, source mRNA, and position in the header. |
| `csv` | The table's columns plus the passenger sequence, with ASCII headers (`Efficacy_Score`, `Delta_G_Score`). |
| `xlsx` | A `shRNA Candidates` sheet with colour scales on the score, GC, and ΔG columns, plus a `Run Info` sheet. |
| `json` | A `meta` object and a `candidates` array with raw numeric fields. |

Every file format carries the same provenance block: tool version, UTC timestamp, input path, a short SHA-256 of the input, and every parameter of the run. It is the `#` comment header in CSV, the `Run Info` sheet in XLSX, and the `meta` object in JSON. `#` lines are ignored on read, so a candidates CSV can be passed back in as a [`--batch`](#a-batch-of-guides-from-a-csv) file.

`-o`/`--output` takes a path without an extension; the extension of the chosen format is appended. Without `-o`, the file is named from the input's stem plus a timestamp and written next to the input, or into the current directory when the input came from stdin. `-o -` writes the data to stdout.

At a terminal, a file format prints the table in addition to writing the file. On a pipe, only the data is written.

`--rna` emits sequences with `U` instead of `T`. Only sequence fields are converted, so backbone names containing a `T` are unaffected.

### Scoring several files at once

Each `INPUT` is loaded, scored, filtered, and exported independently, as if the command had been run once per file.

```bash
# three tables to the terminal
knockknock <mRNAs_A.fasta> <mRNAs_B.fasta> <mRNAs_C.fasta>

# three CSV files to the input dir
knockknock <mRNAs_A.fasta> <mRNAs_B.fasta> <mRNAs_C.fasta> -f csv
```

`-o` is rejected with multiple inputs, because each input produces its own file. `--tui` is rejected, because it displays one dataset. `-` may appear at most once. Multiple inputs always write files instead of streaming to stdout. If two auto-named exports resolve to the same path — same stem, same second — the second is renamed with a `_2` suffix and a warning.

Progress markers and warnings go to stderr, so a redirected data stream carries data only.

## Designing cloning oligos

A backbone added to a scoring run designs oligo pairs for the same candidates that are shown or exported.

```bash
knockknock <mRNAs.fasta> -b pINDUCER
knockknock <mRNAs_A.fasta> <mRNAs_B.fasta> -b pINDUCER -f csv
```

Oligo design is triggered by `-b`/`--backbone` or by any manual insert flag (`--loop`, `--adapter-5`, `--adapter-3`, `--site-5`, `--site-3`). The insert must be defined in one of two ways:

- a named backbone from the database (`--list-backbones`), with `-b`/`--backbone`; or
- all three of `--loop`, `--adapter-5`, and `--adapter-3`. An empty value is valid — `--loop ''` means no loop — but omitting one is an error.

Individual segments of a named backbone can be overridden. `--adapter-5`/`--adapter-3` also accept the names `--flank-5`/`--flank-3`.

```bash
# pINDUCER, but without loop
knockknock <mRNAs.fasta> -b pINDUCER --loop ''

# pINDUCER, but with different ends
knockknock <mRNAs.fasta> -b pINDUCER --site-5 "G/AATTC" --site-3 "GGTAC/C"
```

`--guide-on-5p`/`--guide-on-3p` selects the arm that carries the guide. The default, `--guide-on-3p`, matches every bundled backbone.

`-t`/`--top` has a different default here. It still limits what is shown and exported, and it additionally limits how many candidates get oligos, but its default is every candidate that passed the filters rather than 10. An explicit `-t` applies in either mode.

```bash
# oligos for every passing candidate
knockknock <mRNAs.fasta> -b pINDUCER

# oligos for the top 5 only
knockknock <mRNAs.fasta> -b pINDUCER -t 5
```

### Where the oligos end up

Oligos are added to the requested format, not substituted for it.

| Format | Oligos appear as |
|---|---|
| `table` | A second panel below the candidate table. |
| `fasta` | Two extra records (`fwd`, `rev`) per pair, after the guide records. |
| `json` | An `oligos` array alongside `candidates`. |
| `xlsx` | A second `Oligo Pairs` sheet. |
| `csv` | A separate `*_oligos.csv` file next to the candidates CSV. |

CSV is the exception because one file cannot hold two differently shaped tables. Piping `-f csv` to stdout therefore carries the candidates only; KnockKnock warns and points to `-o` when oligos were designed.

### Defining the backbone ends

`--site-5` and `--site-3` describe the ends of the digested backbone, not the ends of the insert. The oligos are designed to fit them, so a typed sequence reaches an oligo as its reverse complement. Two notations are accepted:

- **cut-site notation**: a recognition sequence with a slash at the top-strand cut position, such as `G/AATTC` (EcoRI) or `CTGCA/G` (PstI). It names the uncut plasmid's top strand and reads the same at either end.
- **a stated overhang**: read 5'→3' along whichever backbone strand protrudes, such as `5'AATT`, `3'TGCA`, or `blunt`. Use this for enzymes that cut outside their recognition sequence (Type IIS, Golden Gate).

A degenerate IUPAC code is accepted in a recognition sequence as long as it does not fall on the overhang. `N/TCGAG` parses, because the overhang is still a definite `TCGA`; `CACNNN/GTG` is rejected, because there the `N` region is the overhang undefined.

### Rebuilding the restriction site

Digestion splits a recognition site between the two fragments. By default the insert carries the overhang only, which leaves a scarred junction that the enzyme can no longer cut. `--keep-sites` adds the insert's half of the recognition sequence, so ligation reforms the site: with pINDUCER's `C/TCGAG` the forward oligo starts `TCGAG…` instead of `TCGA…`, the backbone supplies the leading `C`, and `CTCGAG` is complete. These bases are double-stranded, so both oligos and the annealing Tm grow with them.

`keep_sites` is a stored column of the backbone database, on for `pINDUCER` and off for every other bundled backbone. The flag is three-state: with neither `--keep-sites` nor `--no-keep-sites`, the stored value applies.

```bash
knockknock <mRNAs.fasta> -b pINDUCER                    # on  (stored)
knockknock <mRNAs.fasta> -b pINDUCER --no-keep-sites    # off (this run)
knockknock <mRNAs.fasta> -b pLKO --keep-sites           # on  (this run)
```

Two kinds of junction cannot be rebuilt, and both are reported rather than skipped silently: a junction given as a bare overhang, which has no recognition sequence to restore (normal for Golden Gate backbones such as pHREP, where the Type IIS sites are deliberately not regenerated), and a junction whose degenerate base falls on the insert's share.

A degenerate base is never substituted. The side of the cut it falls on determines the outcome:

| Site and junction | Insert's share | Result |
|---|---|---|
| `CN/TCAGC` at the 5' end | `GC` | Added verbatim. The `N` stays on the vector and is never synthesized. |
| `CN/TCAGC` at the 3' end | `CN` | Nothing added. The site is not padded to length, and the end is reported as unrestorable. |

The warning in the first case states only that the reformation of `CNTCAGC` depends on what the vector carries at the `N`, which KnockKnock cannot know. No base was substituted.

### Mismatch positions

`-m`/`--mismatch` takes one or more 1-based positions, comma- or space-separated, and places a mismatch in the passenger strand at each. The same list is used for the passenger sequence reported by scoring and for the oligo's passenger arm.

Positions count from the base of the stem — the end nearest the flanking sequence, farthest from the loop — regardless of which arm carries the guide. `0` on its own disables mismatch introduction and yields a perfect-complement passenger; `0` combined with any other position is an error. Substitutions are deterministic (`G→T`, `T→G`, `A→C`, `C→A`), so a given guide always yields the same oligo pair.

Without `-m`, a named backbone's stored positions apply instead of the flag's default of `1`. An explicit `-m` overrides both.

### Design warnings

Design warnings are written to stderr, never to the data stream. KnockKnock reports:

- non-directional cloning: both backbone ends present the same overhang, so the insert can ligate in either orientation;
- a backbone enzyme's recognition site occurring inside the insert. This does not affect the ligation, but it matters when re-digesting the finished plasmid;
- any junction that could not be rebuilt under `--keep-sites`.

## Designing oligos from existing guides

`-g`/`--guide` and `--batch` design oligos without scoring a FASTA. The insert is defined exactly as above, by a backbone or by the full manual set of flags. The available formats are narrower, since there is no candidate table: `fasta` (default) and `json`. Neither flag can be combined with `INPUT`, or with the other.

### A single guide

```bash
knockknock -g ATGCATGCATGCATGCATGCA -b pINDUCER
```

The guide must be 19–30 nt. `U` is accepted and read as `T`.

### A batch of guides from a CSV

```bash
knockknock --batch results.csv -b pINDUCER -o oligos
```

The CSV needs a guide-sequence column named `antisense`, `guide`, `guide_sequence`, `Guide_Sequence`, or `Guide Sequence`. The first match wins, checked in that order. Two further columns are read if present:

- a name column (`name`, `Name`, or `shRNA`), used to label each pair in the export. Rows without one become `shRNA-<row number>`.
- a passenger column (`sense`, `passenger`, `passenger_sequence`, or `Passenger Sequence`), currently informational only. The mismatched passenger is always computed, not read from the file.

Lines starting with `#` are ignored, so a candidates CSV from a scoring run can be passed in directly.

Guides in one batch may differ in length. `--mismatch` positions are checked against the shortest guide, so every requested position exists in every assembled oligo.

## Backbone management

`--list-backbones` prints the database: name, both enzymes and ends, both flanks, loop, guide arm, and stored mismatch positions. KnockKnock ships with `pINDUCER`, `pINDUCER_PCR`, `pGIPZ`, `pTRIPZ`, `Open Biosystems`, `pLKO`, and `pHREP`.

Backbones are saved and deleted in the TUI's oligo designer (see [TUI_MANUAL.md](TUI_MANUAL.md#backbone-configuration)). Deletion edits your own copy of the database only, so bundled backbones can be restored:

```bash
knockknock --restore-backbones
```

This adds back any missing bundled backbone and leaves backbones you saved yourself unchanged.

## Information flags

Each of these prints and exits.

- `--version` — the version.
- `--info` — the version, a summary of the bundled reference data (coefficients, ΔG table, dG bins), and the path of the active backbone database.
- `-h`/`--help` — the option summary.

## Option reference

### Input and output

| Option | Effect |
|---|---|
| `INPUT...` | One or more FASTA (or bare-sequence) files to score. `-` reads stdin. Omit to open the TUI at a terminal, or to score piped input. Cannot be combined with `--guide`/`--batch`. |
| `-o`, `--output` | Output path without an extension. `-` writes to stdout. Rejected with multiple inputs. |
| `-f`, `--format` | `table`, `fasta`, `csv`, `xlsx`, `json`. Default: `table` at a terminal, `fasta` when piped. With `--guide`/`--batch`: `fasta` (default) or `json`. |
| `--rna` | Emit sequences as RNA (`U` for `T`). |
| `--tui` | Open the interactive explorer for a single input file. |

### Scoring and filtering

| Option | Effect |
|---|---|
| `-l`, `--length` | Guide length in nt, 20–22 (default 22). Scores are comparable within one length only. |
| `--quality` | Filter preset: `lax`, `standard`, `strict` (default). Individual filter flags override it one at a time. |
| `--gc-min`, `--gc-max` | GC-content window, 0–100 %. |
| `-s`, `--min-efficacy-score` | Reject candidates below this efficacy score. |
| `--remove-motifs` | Reject homopolymeric and low-complexity motifs. |
| `--multitarget` | Require the guide in at least N input mRNAs. Default: all of them. |
| `--sort-dg` | Rank by ΔG score instead of efficacy score. |
| `-t`, `--top` | Keep the top N candidates, and with oligo design on, the N that get oligos. Default: 10 for plain scoring, unlimited once oligo design is triggered. |
| `--seed` | RNG seed. Currently unused and reserved; recorded in export metadata. |

### Oligo design

| Option | Effect |
|---|---|
| `-b`, `--backbone` | Backbone name from the database. Alongside `INPUT`, triggers oligo design for the resulting candidates. |
| `-g`, `--guide` | Design oligos for one 19–30 nt guide, without scoring. |
| `--batch` | Design oligos for a CSV of guides, without scoring. |
| `--loop` | Loop sequence between the two arms. Required, with both adapters, when no backbone is named. |
| `--adapter-5`, `--flank-5` | Constant sequence 5' of the arms. |
| `--adapter-3`, `--flank-3` | Constant sequence 3' of the arms. |
| `--site-5` | Left end of the digested backbone: `G/AATTC`, `5'AATT`, `3'TGCA`, or `blunt`. |
| `--site-3` | Right end of the digested backbone, same notations. |
| `--keep-sites` / `--no-keep-sites` | Carry the insert's half of each recognition sequence, so ligation rebuilds the site. Default: the backbone's stored setting. |
| `-m`, `--mismatch` | Passenger mismatch position(s), 1-based from the base of the stem. `0` alone disables. Unset, a named backbone's stored positions apply. |
| `--guide-on-5p` / `--guide-on-3p` | Which arm carries the guide. Default: `--guide-on-3p`. |

### Information

| Option | Effect |
|---|---|
| `--list-backbones` | Print the backbone database and exit. |
| `--restore-backbones` | Restore any missing bundled backbones and exit. |
| `--info` | Print the version, the bundled reference-data summary, and the active backbone database path, then exit. |
| `--version` | Print the version and exit. |
| `-h`, `--help` | Print the option summary and exit. |
