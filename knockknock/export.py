"""Output formatting: terminal table, CSV, XLSX, JSON."""

import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _distribution_version
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from knockknock.engine import AnalysisParams, Candidate
from knockknock.oligo import OligoPair, fasta_header, maybe_rna, oligo_rows
from knockknock.reference import load_reference_data

console = Console()

#: Label column of the terminal "Filtering Summary" panel, wide enough for the
#: longest label a filter can produce.
_STAT_LABEL_WIDTH = 39


def _tool_version() -> str:
    """The running KnockKnock version, for the provenance block.

    Read from the installed distribution's metadata rather than imported from
    `knockknock.__init__`. That keeps `export` importable even when the package
    root cannot be imported — which is exactly what a stale editable install
    produces, and it should not take down every command just to stamp a version
    string into an export.
    """
    try:
        return _distribution_version("knockknock-shrna")
    except PackageNotFoundError:  # pragma: no cover - uninstalled source tree
        return "0.0.0+unknown"


def resolve_export_path(
    ext: str,
    *,
    user_output: str | None = None,
    reference_path: str | None = None,
    stem: str = "shRNAs",
    label: str | None = None,
) -> Path:
    """Pick an absolute export path — the single source of truth for
    every export site in the CLI and the TUI (main table, oligo designer).

    ``user_output`` is honored as given. Otherwise the file lands beside
    ``reference_path`` (the input the export came from), named after its stem
    plus ``label`` and a timestamp — the same folder the user already has the
    input open in/from, so there's nothing to hunt for afterwards. Without a
    real input file to sit beside (stdin via ``-``, a pasted sequence, or no
    ``reference_path`` at all), the export falls back to the current working
    directory.

    ``label`` marks what kind of export this is (e.g. ``"oligos"``). Two export
    sites can emit the same extension, so without it a candidate table and an
    oligo list exported within the same second would resolve to the same
    filename and one would silently overwrite the other.
    """
    if user_output:
        return Path(user_output).expanduser().with_suffix(ext).resolve()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{label}" if label else ""

    base = stem
    directory = Path.cwd()
    if reference_path and reference_path != "-":
        ref = Path(reference_path).expanduser()
        base = ref.stem or stem
        directory = ref.parent

    return (directory / f"{base}{tag}_{timestamp}{ext}").resolve()


def _metadata(params: AnalysisParams, input_path: str = "") -> dict:
    """Provenance dict embedded in every export.

    Lets a user reproduce a result later: tool version, parameters, input file
    hash, timestamp. Keep all values JSON-serialisable.
    """
    sha = ""
    if input_path and input_path != "-":
        p = Path(input_path)
        if p.is_file():
            sha = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    return {
        "knockknock_version": _tool_version(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input_path": input_path,
        "input_sha256_short": sha,
        "guide_length": params.length,
        "gc_min": params.gc_min,
        "gc_max": params.gc_max,
        "min_efficacy_score": params.min_efficacy_score,
        "remove_motifs": params.remove_motifs,
        "multitarget": params.multitarget,
        "stem_mismatch": list(params.stem_mismatch),
        "sort_by_dg": params.sort_by_dg,
        "rna_output": params.rna_output,
        "seed": params.seed,
    }


def no_candidates_message(stats: dict[str, int], params: AnalysisParams) -> str:
    """Why the candidate list came back empty, phrased for whoever ran it.

    An empty table is indistinguishable from a crash, a mis-parsed file, or a
    tool that simply does nothing — so no output path is allowed to just print
    a header row and stop. The single wording lives here because the CLI
    table, the CLI's export paths and the TUI's picker all have to say the
    same thing.

    Deliberately free of flag names: the TUI reaches this too, and there is no
    ``--quality`` to type in it.
    """
    length = params.length
    if stats.get("windows_scanned", 0) == 0:
        return (
            f"No shRNA candidates: not one {length} nt window fits any input "
            "sequence. Every sequence is shorter than the guide length — give "
            "a longer target, or reduce the guide length."
        )

    unique = stats.get("unique_guides", 0)
    wildcard = stats.get("wildcard_filtered", 0)
    if wildcard and wildcard == unique:
        return (
            f"No shRNA candidates: all {unique} guide windows overlap a "
            "wildcard base. Every defined (A/C/G/T) stretch in the input is "
            f"shorter than the {length} nt guide length, so no guide can be "
            "read out of it."
        )

    reasons = [
        (wildcard, "overlapped a wildcard base"),
        (
            stats.get("gc_filtered", 0),
            f"fell outside GC {params.gc_min:.0f}–{params.gc_max:.0f}%",
        ),
        (
            stats.get("efficacy_score_filtered", 0),
            f"scored below {params.min_efficacy_score}",
        ),
        (stats.get("motif_filtered", 0), "hit an avoided motif"),
        (
            stats.get("target_filtered", 0),
            f"hit fewer than {params.multitarget} of the input mRNAs",
        ),
    ]
    detail = ", ".join(f"{count} {why}" for count, why in reasons if count)
    return (
        f"No shRNA candidates passed the filters: all {unique} unique guide "
        f"windows were rejected ({detail}). Loosen the quality preset, or "
        "widen the GC / efficacy score / multi-target limits."
    )


def format_candidate_cells(
    candidate: Candidate, total_mrnas: int, rna_output: bool
) -> tuple[str, str, str, str, str, str]:
    """The efficacy score/ΔG score/GC/targets cell text shared by every
    candidate *table* view — the CLI's terminal table and the TUI's live
    results table.

    Returns ``(guide_sequence, efficacy_score, delta_g_score, gc_pct, dg,
    targets)``, each already formatted to the precision these two tables
    display. Machine-readable exports (CSV/JSON/XLSX) keep the raw numeric
    fields instead — this is display formatting only, not a serialisation
    format.
    """
    return (
        maybe_rna(candidate.antisense, rna_output),
        f"{candidate.efficacy_score:.4f}",
        f"{candidate.delta_g_score:.1f}",
        f"{candidate.gc_pct:.1f}",
        f"{candidate.dg:.1f}",
        f"{candidate.targets}/{total_mrnas}",
    )


def to_terminal(
    candidates: list[Candidate],
    stats: dict[str, int],
    params: AnalysisParams,
    top_n: int | None = 10,
    total_mrnas: int = 1,
    input_path: str = "",
    oligo_pairs: list[OligoPair] | None = None,
    backbone_name: str | None = None,
) -> None:
    table = Table(
        title=None,
        show_header=True,
        header_style="bold",
        expand=False,
    )
    table.add_column("Rank", justify="right", style="dim")
    table.add_column("Guide Sequence")
    table.add_column("Efficacy Score", justify="right")
    table.add_column("ΔG Score", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("ΔG", justify="right")
    table.add_column("Targets", justify="right")
    table.add_column("Position", justify="right")

    display = candidates[:top_n]

    for rank, c in enumerate(display, 1):
        antisense, efficacy_score, delta_g_score, gc_pct, dg, targets = (
            format_candidate_cells(c, total_mrnas, params.rna_output)
        )
        table.add_row(
            str(rank),
            antisense,
            efficacy_score,
            delta_g_score,
            gc_pct,
            dg,
            targets,
            str(c.start_pos),
        )

    # Summary panel
    summary_lines = [
        f"Input: {input_path}" if input_path else "",
        f"mRNAs: {total_mrnas}",
        f"Candidates: {stats['windows_scanned']} windows "
        f"({stats['unique_guides']} unique), "
        f"{stats['passed']} passed filters",
        f"Guide length: {params.length} nt",
        f"GC range: {params.gc_min:.0f}–{params.gc_max:.0f}%",
        f"Motif filter: {'enabled' if params.remove_motifs else 'disabled'}",
        f"Sort by: {'ΔG score' if params.sort_by_dg else 'efficacy score'}",
    ]
    summary_text = "\n".join(line for line in summary_lines if line)
    console.print(Panel(summary_text, title="Run Summary", border_style="dim"))

    if candidates:
        console.print(table)
    else:
        # Never a bare header row: see `no_candidates_message`.
        console.print(
            Panel(
                no_candidates_message(stats, params),
                title="No shRNA candidates",
                border_style="yellow",
            )
        )

    if top_n is not None and len(candidates) > top_n:
        console.print(
            f"\n[dim]Showing top {top_n} of {len(candidates)} candidates. "
            f"Every export format is capped the same way — raise -t to "
            f"include more, e.g. -t {len(candidates)}.[/dim]"
        )

    # Filtering stats panel. Padded by width rather than by hand-counted
    # spaces: several labels embed a parameter value (efficacy score
    # threshold, multi-target count) whose length isn't known until it's
    # rendered.
    def stat_line(label: str, count: int) -> str:
        return f"{label:<{_STAT_LABEL_WIDTH}}{count:>6}"

    stats_lines = []
    # The wildcard row is only shown when it happened: most inputs carry no
    # wildcards, and a permanent "0" invites the question of what it means.
    if stats.get("wildcard_filtered", 0):
        stats_lines.append(
            stat_line("Rejected by wildcard base:", stats["wildcard_filtered"])
        )
    stats_lines.append(
        stat_line(
            f"Rejected by GC content ({params.gc_min:.0f}–{params.gc_max:.0f}%):",
            stats["gc_filtered"],
        )
    )
    if params.min_efficacy_score is not None:
        stats_lines.append(
            stat_line(
                f"Rejected by minimum efficacy score (< {params.min_efficacy_score}):",
                stats["efficacy_score_filtered"],
            )
        )
    stats_lines.extend(
        [
            stat_line("Rejected by motif filter:", stats["motif_filtered"]),
            stat_line(
                f"Rejected by multi-target (< {params.multitarget}):",
                stats.get("target_filtered", 0),
            ),
            "─" * (_STAT_LABEL_WIDTH + 6),
            stat_line("Passed all filters:", stats["passed"]),
        ]
    )
    console.print(
        Panel("\n".join(stats_lines), title="Filtering Summary", border_style="dim")
    )

    if oligo_pairs:
        oligo_table = Table(
            title=None,
            show_header=True,
            header_style="bold",
            expand=False,
        )
        oligo_table.add_column("Name")
        oligo_table.add_column("Forward (5'->3')")
        oligo_table.add_column("Reverse (5'->3')")
        oligo_table.add_column("Tm", justify="right")
        oligo_table.add_column("Overlap", justify="right")
        for pair in oligo_pairs:
            oligo_table.add_row(
                pair.name,
                maybe_rna(pair.forward, params.rna_output),
                maybe_rna(pair.reverse, params.rna_output),
                f"{pair.tm:.1f}",
                str(pair.overlap),
            )
        title = f"Oligo Pairs ({backbone_name})" if backbone_name else "Oligo Pairs"
        console.print(Panel(oligo_table, title=title, border_style="dim"))


def _csv_text(
    candidates: list[Candidate],
    params: AnalysisParams,
    input_path: str = "",
    total_mrnas: int = 1,
) -> str:
    rows = []
    for rank, c in enumerate(candidates, 1):
        rows.append(
            {
                "Rank": rank,
                "Guide_Sequence": maybe_rna(c.antisense, params.rna_output),
                "Passenger_Sequence": maybe_rna(c.sense, params.rna_output),
                "Efficacy_Score": c.efficacy_score,
                "Delta_G_Score": c.delta_g_score,
                "GC_pct": c.gc_pct,
                "dG": c.dg,
                "Targets": f"{c.targets}/{total_mrnas}",
                "Name": c.name,
                "Start_Pos": c.start_pos,
                "mRNA": c.mrna_id,
            }
        )

    meta = _metadata(params, input_path)
    buf = io.StringIO()
    for k, v in meta.items():
        buf.write(f"# {k}: {v}\n")
    if rows:
        writer = csv.DictWriter(buf, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return buf.getvalue()


def to_csv_text(
    candidates: list[Candidate],
    params: AnalysisParams,
    input_path: str = "",
    total_mrnas: int = 1,
) -> str:
    """CSV rendering as a string.

    Every text format in this module returns a string and lets the caller write
    it; only `to_xlsx` takes a path, because xlsxwriter can only write to one.
    There used to be a second, path-writing `to_csv` as well, with the CLI using
    one and the TUI the other — two ways to do the same thing, each with its own
    encoding and newline handling to get wrong.
    """
    return _csv_text(candidates, params, input_path, total_mrnas)


def _json_payload(
    candidates: list[Candidate],
    params: AnalysisParams,
    input_path: str = "",
    oligo_pairs: list[OligoPair] | None = None,
    backbone_name: str | None = None,
) -> dict:
    rows = []
    for rank, c in enumerate(candidates, 1):
        rows.append(
            {
                "rank": rank,
                "name": c.name,
                "guide_sequence": maybe_rna(c.antisense, params.rna_output),
                "passenger_sequence": maybe_rna(c.sense, params.rna_output),
                "efficacy_score": c.efficacy_score,
                "delta_g_score": c.delta_g_score,
                "gc_pct": c.gc_pct,
                "dg": c.dg,
                "targets": c.targets,
                "start_pos": c.start_pos,
                "mrna_id": c.mrna_id,
            }
        )

    payload = {
        "meta": _metadata(params, input_path),
        "candidates": rows,
    }
    if oligo_pairs is not None:
        payload["oligos"] = oligo_rows(
            oligo_pairs, backbone_name or "Custom", rna=params.rna_output
        )
    return payload


def to_json_text(
    candidates: list[Candidate],
    params: AnalysisParams,
    input_path: str = "",
    oligo_pairs: list[OligoPair] | None = None,
    backbone_name: str | None = None,
) -> str:
    """JSON rendering as a string, for piping to stdout.

    When ``oligo_pairs`` is given, the payload gains a sibling ``"oligos"``
    key alongside ``"candidates"`` — one JSON document, not a second file,
    since JSON (unlike CSV) can hold two distinct sections cleanly.
    """
    payload = _json_payload(candidates, params, input_path, oligo_pairs, backbone_name)
    return json.dumps(payload, indent=2)


def to_fasta(
    candidates: list[Candidate],
    params: AnalysisParams,
    oligo_pairs: list[OligoPair] | None = None,
    backbone_name: str | None = None,
) -> str:
    """FASTA rendering of guide sequences — the default when stdout is piped,
    since that's the format downstream sequence tools (BLAST, primer3,
    aligners) actually expect.

    When ``oligo_pairs`` is given, each pair's fwd/rev oligo is appended as
    two more records after the guide records — FASTA has no concept of a
    second sheet, so distinct records in the same stream are the closest
    match to the other formats' separate oligo section.
    """
    lines = []
    for rank, c in enumerate(candidates, 1):
        antisense = maybe_rna(c.antisense, params.rna_output)
        lines.append(
            f">{c.name} rank={rank} efficacy_score={c.efficacy_score:.4f} "
            f"gc={c.gc_pct:.1f} mRNA={c.mrna_id} pos={c.start_pos}"
        )
        lines.append(antisense)
    if oligo_pairs:
        name = backbone_name or "Custom"
        for pair in oligo_pairs:
            lines.append(fasta_header(pair.name, "fwd", name))
            lines.append(maybe_rna(pair.forward, params.rna_output))
            lines.append(fasta_header(pair.name, "rev", name))
            lines.append(maybe_rna(pair.reverse, params.rna_output))
    return "\n".join(lines) + ("\n" if lines else "")


def _apply_scale(
    worksheet,
    col: int,
    n_rows: int,
    *,
    min_value: float,
    mid_value: float,
    max_value: float,
    min_color: str | None = None,
    mid_color: str | None = None,
    max_color: str | None = None,
) -> None:
    """Apply a 3-color conditional-format scale to one data column.

    Shared by every `to_xlsx` column scale below, which otherwise differ
    only in column index and these min/mid/max value/color numbers.
    """
    options = {
        "type": "3_color_scale",
        "min_value": min_value,
        "mid_value": mid_value,
        "max_value": max_value,
        "min_type": "num",
        "mid_type": "num",
        "max_type": "num",
    }
    if min_color is not None:
        options["min_color"] = min_color
    if mid_color is not None:
        options["mid_color"] = mid_color
    if max_color is not None:
        options["max_color"] = max_color
    worksheet.conditional_format(1, col, n_rows, col, options)


def to_xlsx(
    candidates: list[Candidate],
    params: AnalysisParams,
    output_path: Path,
    total_mrnas: int = 1,
    input_path: str = "",
    oligo_pairs: list[OligoPair] | None = None,
    backbone_name: str | None = None,
) -> None:
    rows = []
    for rank, c in enumerate(candidates, 1):
        rows.append(
            {
                "Rank": rank,
                "Guide": maybe_rna(c.antisense, params.rna_output),
                "Passenger": maybe_rna(c.sense, params.rna_output),
                "Efficacy Score": c.efficacy_score,
                "ΔG Score": c.delta_g_score,
                "GC%": c.gc_pct,
                "ΔG": c.dg,
                "Targets": f"{c.targets}/{total_mrnas}",
                "Position": c.start_pos,
                "mRNA": c.mrna_id,
            }
        )
    columns = [
        "Rank",
        "Guide",
        "Passenger",
        "Efficacy Score",
        "ΔG Score",
        "GC%",
        "ΔG",
        "Targets",
        "Position",
        "mRNA",
    ]
    # Write a usable workbook even when no candidate passes the filters.
    # An empty sheet lets users distinguish "no hits" from a failed export.
    df = pd.DataFrame(rows, columns=columns)
    excel_writer = pd.ExcelWriter(str(output_path), engine="xlsxwriter")
    sheet_name = "shRNA Candidates"
    df.to_excel(excel_writer, sheet_name=sheet_name, index=False)

    worksheet = excel_writer.sheets[sheet_name]

    for idx, col in enumerate(df.columns):
        series = df[col]
        max_len = max([len(str(series.name))] + [len(str(value)) for value in series]) + 1
        worksheet.set_column(idx, idx, max_len)

    # Efficacy score / ΔG score colour scales are derived from the reference
    # data so they stay in sync if the coefficient tables are ever updated.
    ref = load_reference_data()
    efficacy_score_min = float(ref.coeffs[["A", "C", "G", "T"]].min(axis=1).sum())
    efficacy_score_max = float(ref.coeffs[["A", "C", "G", "T"]].max(axis=1).sum())
    delta_g_score_min = float(ref.dG_bins.min(axis=1).sum())
    delta_g_score_max = float(ref.dG_bins.max(axis=1).sum())
    delta_g_score_mid = (delta_g_score_min + delta_g_score_max) / 2

    n_rows = len(candidates)
    if n_rows:
        _apply_scale(
            worksheet,
            3,
            n_rows,
            min_value=efficacy_score_min,
            mid_value=0,
            max_value=efficacy_score_max,
        )
        _apply_scale(
            worksheet,
            4,
            n_rows,
            min_value=delta_g_score_min,
            mid_value=delta_g_score_mid,
            max_value=delta_g_score_max,
        )
        _apply_scale(
            worksheet,
            5,
            n_rows,
            min_value=0,
            mid_value=50,
            max_value=100,
            min_color="#63BE7B",
            max_color="#F8696B",
        )
        _apply_scale(
            worksheet,
            6,
            n_rows,
            min_value=-40,
            mid_value=-33,
            max_value=-20,
            min_color="#F8696B",
            mid_color="#63BE7B",
            max_color="#F8696B",
        )

    # Oligo Pairs sheet: a second, distinct sheet rather than more columns on
    # the candidates sheet — keeps the two concerns (scoring vs. cloning
    # oligos) visually and structurally separate, same as the terminal panel.
    if oligo_pairs:
        oligo_columns = [
            "name",
            "backbone",
            "forward",
            "reverse",
            "guide",
            "passenger",
            "overlap",
            "tm",
        ]
        oligo_df = pd.DataFrame(
            oligo_rows(oligo_pairs, backbone_name or "Custom", rna=params.rna_output),
            columns=oligo_columns,
        )
        oligo_sheet = "Oligo Pairs"
        oligo_df.to_excel(excel_writer, sheet_name=oligo_sheet, index=False)
        oligo_ws = excel_writer.sheets[oligo_sheet]
        for idx, col in enumerate(oligo_df.columns):
            series = oligo_df[col]
            max_len = (
                max([len(str(series.name))] + [len(str(value)) for value in series]) + 1
            )
            oligo_ws.set_column(idx, idx, max_len)

    # Provenance: key/value sheet so the workbook is self-describing.
    meta_df = pd.DataFrame(
        list(_metadata(params, input_path).items()),
        columns=["key", "value"],
    )
    meta_df.to_excel(excel_writer, sheet_name="Run Info", index=False)
    meta_ws = excel_writer.sheets["Run Info"]
    meta_ws.set_column(0, 0, max([len(str(v)) for v in meta_df["key"]] + [5]) + 1)
    meta_ws.set_column(1, 1, max([len(str(v)) for v in meta_df["value"]] + [5]) + 1)

    excel_writer.close()


def write_candidates_export(
    fmt: str,
    candidates: list[Candidate],
    params: AnalysisParams,
    out_path: Path,
    *,
    input_path: str = "",
    total_mrnas: int = 1,
) -> None:
    """Write a plain candidate-table export (no oligo section) to ``out_path``.

    The single place a format string (``"csv"``/``"xlsx"``/``"json"``/
    ``"fasta"``) turns into a call to the right `to_*` function above —
    `tui.py`'s explorer screen and `oligo_tui.py`'s oligo designer both
    dispatch their candidate export through this, so a screen exporting
    candidates can't drift from what the other one writes for the same
    format.
    """
    if fmt == "xlsx":
        to_xlsx(
            candidates, params, out_path, total_mrnas=total_mrnas, input_path=input_path
        )
        return
    if fmt == "csv":
        text = to_csv_text(candidates, params, input_path=input_path, total_mrnas=total_mrnas)
    elif fmt == "json":
        text = to_json_text(candidates, params, input_path=input_path)
    else:  # fasta
        text = to_fasta(candidates, params)
    out_path.write_text(text, encoding="utf-8")
