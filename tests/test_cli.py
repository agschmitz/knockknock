"""Integration tests for the Typer CLI."""

import json
import re
import sys
from pathlib import Path

import pytest
import typer
from Bio.Seq import Seq
from typer.testing import CliRunner

from knockknock import __version__
from knockknock.backbone import Backbone, save_backbone
from knockknock.cli import app
from knockknock.cutsite import parse_overhang

FIXTURE = Path(__file__).parent / "data" / "sample.fa"
runner = CliRunner()


def _stdout_wrapper_class() -> type:
    """The concrete class CliRunner wraps `sys.stdout` in during `invoke()`.

    Resolved dynamically via a throwaway probe command rather than importing
    it by name (e.g. `click.testing._NamedTextIOWrapper`) — that private
    class's import path isn't stable across click/typer versions, and in at
    least one real environment `click` wasn't independently importable even
    though `typer.testing.CliRunner` (which depends on it) worked fine.
    Reading `type(sys.stdout)` from inside an actual invocation sidesteps
    needing to name the module at all.
    """
    probe = typer.Typer()
    captured: list[type] = []

    @probe.command()
    def _probe() -> None:
        captured.append(type(sys.stdout))

    CliRunner().invoke(probe, [])
    return captured[0]


_STDOUT_WRAPPER_CLASS = _stdout_wrapper_class()


def _passed(stdout: str) -> int:
    m = re.search(r"Passed all filters:\s+(\d+)", stdout)
    assert m, f"could not find 'Passed all filters' in:\n{stdout}"
    return int(m.group(1))


# ---------------------------------------------------------------------------
# `run` command
# ---------------------------------------------------------------------------


def test_run_table_format_smoke():
    """Explicit `-f table` exits 0 and prints the rich summary + table."""
    result = runner.invoke(app, ["run", str(FIXTURE), "-t", "5", "-f", "table"])
    assert result.exit_code == 0
    assert "Run Summary" in result.stdout
    assert "Passed all filters" in result.stdout


def test_run_csv_writes_file(tmp_path: Path):
    """`-f csv -o <path>` writes a `.csv` file (extension added) with the
    provenance comment header and a Guide_Sequence column."""
    out = tmp_path / "out"
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-f", "csv", "-o", str(out), "--multitarget", "1"]
    )
    assert result.exit_code == 0
    csv_path = out.with_suffix(".csv")
    assert csv_path.exists()
    content = csv_path.read_text()
    assert "Guide_Sequence" in content
    assert content.startswith("# knockknock_version:")


def test_run_json_writes_file_with_meta(tmp_path: Path):
    """`-f json -o <path>` writes a `.json` file whose top level has both
    `meta` (with the resolved input_path) and a non-empty `candidates` list."""
    out = tmp_path / "out"
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-f", "json", "-o", str(out), "--multitarget", "1"]
    )
    assert result.exit_code == 0
    json_path = out.with_suffix(".json")
    data = json.loads(json_path.read_text())
    assert "meta" in data and "candidates" in data
    assert data["meta"]["input_path"] == str(FIXTURE)
    assert len(data["candidates"]) > 0


def test_run_default_format_is_fasta_when_piped():
    """CliRunner's captured stdout is not a tty, so with no -f the tool
    should default to FASTA (piped) rather than the rich table."""
    result = runner.invoke(app, ["run", str(FIXTURE), "-t", "5", "--multitarget", "1"])
    assert result.exit_code == 0
    assert result.stdout.startswith(">")
    assert "Run Summary" not in result.stdout


def test_run_default_format_is_table_at_terminal(monkeypatch):
    """Counterpart to the piped-default test: with stdout patched to report
    as a tty and no explicit -f, the table format is picked, not FASTA."""
    # CliRunner wraps stdout in an in-memory buffer that always reports
    # isatty()==False; patch its class (see `_stdout_wrapper_class` above)
    # to simulate an interactive terminal.
    monkeypatch.setattr(_STDOUT_WRAPPER_CLASS, "isatty", lambda self: True)
    result = runner.invoke(app, ["run", str(FIXTURE), "-t", "5"])
    assert result.exit_code == 0
    assert "Run Summary" in result.stdout


def test_run_explicit_table_format_while_piped():
    """An explicit -f table should win even though stdout isn't a tty."""
    result = runner.invoke(app, ["run", str(FIXTURE), "-t", "5", "-f", "table"])
    assert result.exit_code == 0
    assert "Run Summary" in result.stdout


def test_run_fasta_to_stdout_content():
    """FASTA output is a strict header/sequence pair per candidate, with no
    extra table/summary text mixed in."""
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "5", "-f", "fasta", "--multitarget", "1"]
    )
    assert result.exit_code == 0
    assert result.stdout.startswith(">")
    lines = result.stdout.splitlines()
    assert lines[1] and not lines[1].startswith(">")
    assert len(lines) == 10  # 5 records x (header + sequence)


def test_run_fasta_respects_top_n(tmp_path: Path):
    """-t caps the number of FASTA records emitted, same as the table."""
    unbounded = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "1000", "-f", "fasta", "--multitarget", "1"]
    )
    capped = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "3", "-f", "fasta", "--multitarget", "1"]
    )
    assert unbounded.exit_code == 0 and capped.exit_code == 0
    n_unbounded = sum(1 for ln in unbounded.stdout.splitlines() if ln.startswith(">"))
    n_capped = sum(1 for ln in capped.stdout.splitlines() if ln.startswith(">"))
    assert n_capped == 3
    assert n_unbounded > n_capped


def test_run_csv_respects_top_n():
    """Every export format is capped by -t the same way the table is —
    copy-pasted/exported candidates always match what -t asked for."""
    capped = runner.invoke(
        app,
        ["run", str(FIXTURE), "-t", "1", "-f", "csv", "-o", "-", "--multitarget", "1"],
    )
    uncapped = runner.invoke(
        app,
        [
            "run",
            str(FIXTURE),
            "-t",
            "1000",
            "-f",
            "csv",
            "-o",
            "-",
            "--multitarget",
            "1",
        ],
    )
    assert capped.exit_code == 0 and uncapped.exit_code == 0

    def _data_rows(stdout: str) -> list[str]:
        return [
            ln
            for ln in stdout.splitlines()
            if ln and not ln.startswith("#") and not ln.startswith("Rank,")
        ]

    assert len(_data_rows(capped.stdout)) == 1
    assert len(_data_rows(uncapped.stdout)) > 1


def test_run_csv_dash_forces_stdout_even_at_terminal(monkeypatch):
    """`-o -` always means stdout, even when isatty() is patched to True —
    an explicit dash outranks the terminal-detection default."""
    monkeypatch.setattr(_STDOUT_WRAPPER_CLASS, "isatty", lambda self: True)
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "5", "-f", "csv", "-o", "-", "--multitarget", "1"]
    )
    assert result.exit_code == 0
    assert "Guide_Sequence" in result.stdout
    assert "Run Summary" not in result.stdout


def test_run_xlsx_writes_file(tmp_path: Path):
    """`-f xlsx -o <path>` writes a non-empty `.xlsx` file."""
    out = tmp_path / "out"
    result = runner.invoke(app, ["run", str(FIXTURE), "-f", "xlsx", "-o", str(out)])
    assert result.exit_code == 0
    xlsx_path = out.with_suffix(".xlsx")
    assert xlsx_path.exists() and xlsx_path.stat().st_size > 0


def test_run_quality_strict_smaller_than_lax():
    """The `strict` quality preset must never pass more candidates than `lax`."""
    res_lax = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "0", "-f", "table", "--quality", "lax"]
    )
    res_strict = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "0", "-f", "table", "--quality", "strict"]
    )
    assert res_lax.exit_code == 0 and res_strict.exit_code == 0
    assert _passed(res_strict.stdout) <= _passed(res_lax.stdout)


def test_run_quality_explicit_flag_overrides_preset():
    """--gc-min 0 should override the strict preset's 30."""
    base = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "0", "-f", "table", "--quality", "strict"]
    )
    override = runner.invoke(
        app,
        [
            "run",
            str(FIXTURE),
            "-t",
            "0",
            "-f",
            "table",
            "--quality",
            "strict",
            "--gc-min",
            "0",
            "--gc-max",
            "100",
        ],
    )
    assert base.exit_code == 0 and override.exit_code == 0
    assert _passed(override.stdout) >= _passed(base.stdout)


def test_run_invalid_fasta_alphabet(tmp_path: Path):
    """Non-nucleotide characters in the FASTA sequence exit 1 with an
    'invalid characters' message rather than crashing downstream."""
    bad = tmp_path / "bad.fa"
    bad.write_text(">x\nATGCXY\n")
    result = runner.invoke(app, ["run", str(bad), "-t", "0"])
    assert result.exit_code == 1
    assert "invalid characters" in result.stdout


def test_run_invalid_quality_value():
    """An unrecognized --quality name exits 1 with an 'Unknown --quality' error."""
    result = runner.invoke(app, ["run", str(FIXTURE), "-t", "0", "--quality", "bogus"])
    assert result.exit_code == 1
    assert "Unknown --quality" in result.stdout


def test_run_lowercase_fasta_accepted(tmp_path: Path):
    """Lowercase sequence input is valid DNA, not rejected as invalid alphabet."""
    lo = tmp_path / "lo.fa"
    lo.write_text(">x\natgcgtagctagctagctagctgat\n")
    result = runner.invoke(app, ["run", str(lo), "-t", "0"])
    assert result.exit_code == 0


@pytest.mark.parametrize("sequence", ["ATGCN" * 5, "AUGCAUGCAUGCAUGCAUGCAUGC"])
def test_run_accepts_wildcards_and_rna_without_crashing(tmp_path: Path, sequence: str):
    """IUPAC wildcard bases and RNA (U-containing) input both run to
    completion instead of being rejected up front."""
    fasta = tmp_path / "input.fa"
    fasta.write_text(f">x\n{sequence}\n")
    result = runner.invoke(app, ["run", str(fasta), "-t", "0", "-f", "table"])
    assert result.exit_code == 0


def test_run_warns_about_wildcards_and_keeps_the_defined_stretches(tmp_path: Path):
    """An N run must cost only the guides that overlap it — with a warning
    saying so, rather than the whole record being rejected."""
    fasta = tmp_path / "input.fa"
    defined = "ATGCGTAGCTAGCTAGCTAGCTGAT"
    fasta.write_text(f">x\n{defined}NNNNN{defined[::-1]}\n")

    result = runner.invoke(app, ["run", str(fasta), "-t", "5", "-f", "table"])

    assert result.exit_code == 0
    assert "wildcard base(s)" in result.output
    assert "Rejected by wildcard base" in result.output
    assert _passed(result.output) > 0


def test_run_reports_no_candidates_instead_of_an_empty_table(tmp_path: Path):
    """Every defined stretch shorter than the guide length: nothing to design."""
    fasta = tmp_path / "input.fa"
    fasta.write_text(">x\n" + "ATGCGTAGCN" * 5 + "\n")

    result = runner.invoke(app, ["run", str(fasta), "-f", "table"])

    assert result.exit_code == 0
    assert "No shRNA candidates" in result.output
    assert "wildcard" in result.output


def test_run_reports_no_candidates_when_filters_reject_everything():
    """An impossibly high --min-efficacy-score still exits 0, reporting the
    filters-rejected-everything cause rather than an empty table."""
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-f", "table", "--min-efficacy-score", "9999"]
    )

    assert result.exit_code == 0
    assert "No shRNA candidates passed the filters" in result.output


def test_run_reports_no_candidates_on_a_piped_export(tmp_path: Path):
    """An empty FASTA on stdout is otherwise completely silent."""
    fasta = tmp_path / "input.fa"
    fasta.write_text(">x\n" + "ATGCGTAGCN" * 5 + "\n")

    result = runner.invoke(app, ["run", str(fasta), "-f", "fasta", "-o", "-"])

    assert result.exit_code == 0
    assert "No shRNA candidates" in result.output


@pytest.mark.parametrize(
    "args, message",
    [
        (["--length", "23"], "--length must be between 20 and 22"),
        (["--length", "19"], "--length must be between 20 and 22"),
        (["--top", "-1"], "--top must be zero or greater"),
        (["--multitarget", "0"], "--multitarget must be at least 1"),
        (["--mismatch", "99"], "--mismatch position(s) [99] are outside the 1–22 guide range"),
        (["--mismatch", "0,3"], "--mismatch position 0 (no mismatch) cannot be combined"),
    ],
)
def test_run_rejects_invalid_numeric_options(args, message):
    """Out-of-range numeric options (--length, --top, --multitarget,
    --mismatch) each exit 1 with their specific validation message."""
    result = runner.invoke(app, ["run", str(FIXTURE), *args])
    assert result.exit_code == 1
    assert message in result.stdout


def test_run_rejects_duplicate_fasta_ids(tmp_path: Path):
    """Two FASTA records sharing the same identifier exit 1 rather than
    silently overwriting one another downstream."""
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">same\nATGCATGCATGCATGCATGCATGC\n>same\nATGCATGCATGCATGCATGCATGC\n")
    result = runner.invoke(app, ["run", str(fasta), "-t", "0"])
    assert result.exit_code == 1
    assert "Duplicate FASTA identifier" in result.stdout


def test_run_rejects_unknown_format():
    """An unrecognized -f value exits 1 with an 'Unsupported format' error."""
    result = runner.invoke(app, ["run", str(FIXTURE), "-f", "bogus"])
    assert result.exit_code == 1
    assert "Unsupported format" in result.stdout


def test_run_rejects_gc_min_gte_gc_max():
    """--gc-min at or above --gc-max exits 1 with a 'GC limits must satisfy' error."""
    result = runner.invoke(app, ["run", str(FIXTURE), "--gc-min", "80", "--gc-max", "20"])
    assert result.exit_code == 1
    assert "GC limits must satisfy" in result.stdout


def test_run_multitarget_clamped_to_input_count():
    """A single-mRNA input with --multitarget 5 should warn and clamp rather
    than reject every candidate outright."""
    result = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "0", "-f", "table", "--multitarget", "5"]
    )
    assert result.exit_code == 0
    assert "multitarget reduced to" in result.stdout


def test_run_default_multitarget_requires_every_input_mrna():
    """With no --multitarget, a guide must be shared by every input mRNA —
    same default as the TUI's "target all inputs" checkbox. FIXTURE's two
    records share no common guide, so the default run reports zero
    candidates while an explicit --multitarget 1 finds some."""
    default_run = runner.invoke(app, ["run", str(FIXTURE), "-t", "0", "-f", "table"])
    assert default_run.exit_code == 0
    assert _passed(default_run.stdout) == 0

    single_target = runner.invoke(
        app, ["run", str(FIXTURE), "-t", "0", "-f", "table", "--multitarget", "1"]
    )
    assert single_target.exit_code == 0
    assert _passed(single_target.stdout) > 0


# ---------------------------------------------------------------------------
# `run` with multiple INPUT files (batch mode)
# ---------------------------------------------------------------------------


def _write_gene_fastas(tmp_path: Path) -> tuple[Path, Path]:
    gene_a = tmp_path / "geneA.fa"
    gene_a.write_text(">geneA\nATGCGTAGCTAGCTAGCTAGCTGATCGTAGCTAGCTAGCTAGCTGGTACGTAGCT\n")
    gene_b = tmp_path / "geneB.fa"
    gene_b.write_text(
        ">geneB\nGCATCGATCGTAGCTAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT\n"
    )
    return gene_a, gene_b


def test_run_multiple_inputs_writes_one_file_per_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Two INPUT files with no -o each get their own auto-named CSV, scored
    independently — this is the batch-design use case: pass every gene's
    FASTA in one call and get one export per gene."""
    gene_a, gene_b = _write_gene_fastas(tmp_path)
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(
        app,
        ["run", str(gene_a), str(gene_b), "-f", "csv", "--multitarget", "1"],
    )
    assert result.exit_code == 0
    assert "geneA.fa" in result.output and "geneB.fa" in result.output

    csvs = sorted(tmp_path.glob("*.csv"))
    assert len(csvs) == 2
    stems = {p.stem.split("_")[0] for p in csvs}
    assert stems == {"geneA", "geneB"}
    for csv_path in csvs:
        content = csv_path.read_text()
        assert "Guide_Sequence" in content
        assert content.startswith("# knockknock_version:")


def test_run_multiple_inputs_rejects_explicit_output(tmp_path: Path):
    """-o names a single path, which can't serve two independent exports —
    multiple INPUTs plus an explicit -o must exit 1, not silently let the
    second file overwrite the first."""
    gene_a, gene_b = _write_gene_fastas(tmp_path)
    result = runner.invoke(
        app,
        ["run", str(gene_a), str(gene_b), "-o", str(tmp_path / "out"), "-f", "csv"],
    )
    assert result.exit_code == 1
    assert "cannot be used with multiple input files" in result.stdout


def test_run_multiple_inputs_rejects_tui(tmp_path: Path):
    """--tui only makes sense for one candidate set at a time."""
    gene_a, gene_b = _write_gene_fastas(tmp_path)
    result = runner.invoke(app, ["run", str(gene_a), str(gene_b), "--tui"])
    assert result.exit_code == 1
    assert "--tui can only be used with a single input file" in result.stdout


def test_run_multiple_inputs_dash_output_forces_stdout(tmp_path: Path):
    """`-o -` is not a real path, so it's exempt from the multi-input -o ban
    and instead forces every file's export to stdout, one after another."""
    gene_a, gene_b = _write_gene_fastas(tmp_path)
    result = runner.invoke(
        app,
        [
            "run",
            str(gene_a),
            str(gene_b),
            "-o",
            "-",
            "-f",
            "fasta",
            "--multitarget",
            "1",
        ],
    )
    assert result.exit_code == 0
    assert result.stdout.count(">geneA") >= 1
    assert result.stdout.count(">geneB") >= 1
    assert not list(tmp_path.glob("*.fasta"))


def test_run_rejects_stdin_given_twice():
    """'-' can only be read once; a second '-' among multiple INPUTs exits 1
    instead of trying to re-read an already-drained stdin."""
    result = runner.invoke(app, ["run", "-", "-"])
    assert result.exit_code == 1
    assert "can only be given once" in result.stdout


def test_run_multiple_inputs_disambiguates_same_stem_export(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Two different input files in the same folder that share a filename
    stem (e.g. differing only by extension) must not silently overwrite each
    other's auto-named export, which now lands beside each input — the
    second gets a disambiguating suffix instead."""
    monkeypatch.chdir(tmp_path)
    gene_a = tmp_path / "gene.fa"
    gene_a.write_text(">geneA\nATGCGTAGCTAGCTAGCTAGCTGATCGTAGCTAGCTAGCTAGCTGGTACGTAGCT\n")
    gene_b = tmp_path / "gene.fna"
    gene_b.write_text(
        ">geneB\nGCATCGATCGTAGCTAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT\n"
    )

    result = runner.invoke(
        app,
        ["run", str(gene_a), str(gene_b), "-f", "fasta", "--multitarget", "1"],
    )
    assert result.exit_code == 0
    assert "collided with an earlier file" in result.output

    fastas = sorted(tmp_path.glob("gene*.fasta"))
    assert len(fastas) == 2
    assert fastas[1].stem.endswith("_2")


def test_run_empty_fasta_file(tmp_path: Path):
    """A zero-byte input file exits 1 with 'No sequences found'."""
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    result = runner.invoke(app, ["run", str(empty)])
    assert result.exit_code == 1
    assert "No sequences found" in result.stdout


def test_run_empty_sequence_record(tmp_path: Path):
    """A FASTA record with a header but no sequence exits 1, naming that
    record as empty, even when other records in the file are valid."""
    fasta = tmp_path / "empty_record.fa"
    fasta.write_text(">x\n\n>y\nATGCATGCATGCATGCATGCATGC\n")
    result = runner.invoke(app, ["run", str(fasta)])
    assert result.exit_code == 1
    assert "is empty" in result.stdout


def test_run_stdin_fasta():
    """`INPUT` of `-` reads a FASTA-formatted sequence from stdin."""
    result = runner.invoke(
        app, ["run", "-", "-t", "0"], input=">x\nATGCATGCATGCATGCATGCATGC\n"
    )
    assert result.exit_code == 0


def test_run_stdin_raw_sequence():
    """A pasted sequence with no '>' header is accepted, not just FASTA."""
    result = runner.invoke(
        app, ["run", "-", "-t", "0"], input="ATGCATGCATGCATGCATGCATGC\n"
    )
    assert result.exit_code == 0


def test_run_no_input_argument_reads_stdin():
    """Omitting INPUT entirely is equivalent to passing '-'."""
    result = runner.invoke(
        app, ["run", "-t", "0"], input=">x\nATGCATGCATGCATGCATGCATGC\n"
    )
    assert result.exit_code == 0


def test_info_flag_smoke():
    """`--info` exits 0 and reports the reference data provenance (coeffs.csv,
    Matveeva, Xia), plus the version."""
    result = runner.invoke(app, ["--info"])
    assert result.exit_code == 0
    assert "coeffs.csv" in result.stdout
    assert "Matveeva" in result.stdout
    assert "Xia" in result.stdout
    assert __version__ in result.stdout


def test_version_flag_prints_only_the_version():
    """`--version` exits 0 and prints just the version, not the reference
    data summary `--info` shows."""
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "coeffs.csv" not in result.stdout


def test_dash_h_is_a_help_alias():
    """`-h` behaves like `--help`, both at the top level and on `run`."""
    top = runner.invoke(app, ["-h"])
    long_form = runner.invoke(app, ["--help"])
    assert top.exit_code == 0
    assert top.stdout == long_form.stdout

    run_short = runner.invoke(app, ["run", "-h"])
    assert run_short.exit_code == 0
    assert "Usage:" in run_short.stdout


# ---------------------------------------------------------------------------
# Standalone oligo design (--guide/--batch/--list-backbones, no INPUT)
# ---------------------------------------------------------------------------


def test_oligo_list_backbones():
    """`--list-backbones` exits 0 and prints the known backbone presets."""
    result = runner.invoke(app, ["--list-backbones"])
    assert result.exit_code == 0
    assert "pINDUCER" in result.stdout


def test_restore_backbones_reports_nothing_missing_by_default():
    """`--restore-backbones` on an untouched database says so and exits 0."""
    result = runner.invoke(app, ["--restore-backbones"])
    assert result.exit_code == 0
    assert "Nothing to restore" in result.stdout


def test_restore_backbones_brings_back_a_deleted_bundled_backbone():
    """Deleting a bundled backbone, then `--restore-backbones`, brings it back —
    this is the flow the TUI's delete-confirmation dialog points users at."""
    from knockknock.backbone import delete_backbone

    delete_backbone("pINDUCER")
    result = runner.invoke(app, ["--restore-backbones"])
    assert result.exit_code == 0
    assert "pINDUCER" in result.stdout

    listed = runner.invoke(app, ["--list-backbones"])
    assert "pINDUCER" in listed.stdout


def test_oligo_single_guide_with_backbone():
    """A single -g guide with -b backbone exits 0 and prints both the
    fwd and rev oligo strands."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
        ],
    )
    assert result.exit_code == 0
    assert "_fwd_" in result.stdout and "_rev_" in result.stdout


@pytest.mark.parametrize("flag", ["--loop", "--adapter-5", "--adapter-3"])
def test_oligo_rejects_non_dna_override(flag):
    """A non-DNA character in a flank/loop override must be reported, not crash."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            flag,
            "ACG1",
        ],
    )
    assert result.exit_code == 1
    assert "non-DNA character" in result.stdout


@pytest.mark.parametrize("flag", ["--loop", "--adapter-5", "--adapter-3"])
def test_oligo_rejects_degenerate_override(flag):
    """An IUPAC degenerate code (e.g. N) in a flank/loop override exits 1
    with a 'degenerate IUPAC code' error, distinct from a non-DNA character."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            flag,
            "ACGN",
        ],
    )
    assert result.exit_code == 1
    assert "degenerate IUPAC code" in result.stdout


def test_oligo_batch_csv(tmp_path: Path):
    """`--batch <csv> -o <path>` writes one fwd/rev oligo pair per
    row to a `.fasta` file, headered by each row's name column."""
    csv_path = tmp_path / "batch.csv"
    csv_path.write_text(
        "name,antisense\nsh-1,ATGCATGCATGCATGCATGCAT\nsh-2,GCTAGCTAGCTAGCTAGCTAGC\n"
    )
    out = tmp_path / "oligos"
    result = runner.invoke(
        app,
        [
            "--batch",
            str(csv_path),
            "-b",
            "pINDUCER",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    out_path = out.with_suffix(".fasta")
    assert out_path.exists()
    text = out_path.read_text()
    assert ">sh-1_fwd" in text and ">sh-2_rev" in text


def test_oligo_accepts_candidate_csv_export(tmp_path: Path):
    """A CSV exported by scoring a FASTA can be fed straight back in as
    `--batch`, without reshaping columns first."""
    candidates = tmp_path / "candidates"
    exported = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-f",
            "csv",
            "-o",
            str(candidates),
            "--multitarget",
            "1",
        ],
    )
    assert exported.exit_code == 0

    result = runner.invoke(
        app,
        [
            "--batch",
            str(candidates.with_suffix(".csv")),
            "-b",
            "pINDUCER",
        ],
    )
    assert result.exit_code == 0
    assert "_fwd_" in result.stdout


@pytest.mark.parametrize("guide", ["ATGCN" * 4 + "AT", "ATGC"])
def test_oligo_rejects_invalid_guide(guide: str):
    """A guide that's too short or wildcard-laden exits 1 with a 'Guide' error."""
    result = runner.invoke(app, ["--guide", guide, "-b", "pINDUCER"])
    assert result.exit_code == 1
    assert "Guide" in result.stdout


# A backbone (or an explicit loop + both flanks) is mandatory. Without one there
# is no design strategy, and the tool used to silently emit a loop-less blunt
# guide-fused-to-passenger oligo that looked ready to order.


@pytest.mark.parametrize(
    "args",
    [
        pytest.param(["-g", "ATGCATGCATGCATGCATGCAT"], id="single"),
        pytest.param(["--batch", "BATCH_CSV"], id="batch"),
    ],
)
def test_oligo_requires_a_design_strategy(tmp_path: Path, args):
    """Neither a single guide nor a batch run without a backbone (or a full
    loop+adapter override): both exit 1 naming --backbone and --loop as
    the ways to fix it, so the guide-fused-to-passenger regression can't
    silently ship again."""
    csv_path = tmp_path / "batch.csv"
    csv_path.write_text("name,antisense\nsh-1,ATGCATGCATGCATGCATGCAT\n")
    args = [str(csv_path) if a == "BATCH_CSV" else a for a in args]

    result = runner.invoke(app, args)

    assert result.exit_code == 1
    assert "No oligo design strategy" in result.stdout
    # The message must name what to do next, not just what went wrong.
    assert "--backbone" in result.stdout
    assert "--loop" in result.stdout


def test_oligo_partial_strategy_names_only_the_missing_flags():
    """Supplying only --loop (no backbone) exits 1 naming just the still-missing
    --adapter-5/--adapter-3 flags, not --loop itself."""
    result = runner.invoke(
        app,
        ["-g", "ATGCATGCATGCATGCATGCAT", "--loop", "TAGTGAAGCCACAGATGTA"],
    )
    assert result.exit_code == 1
    assert "--adapter-5" in result.stdout and "--adapter-3" in result.stdout


def test_oligo_accepts_explicit_empty_segments_as_a_strategy():
    """An empty value is a deliberate "this segment is absent", not a missing
    answer — so it satisfies the requirement rather than tripping it."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "--loop",
            "TAGTGAAGCCACAGATGTA",
            "--adapter-5",
            "",
            "--adapter-3",
            "",
        ],
    )
    assert result.exit_code == 0
    assert "_fwd_" in result.stdout and "_rev_" in result.stdout


def test_oligo_empty_override_clears_a_backbone_segment():
    """`--loop ''` must override the backbone's loop, not be ignored as falsy."""
    with_loop = runner.invoke(
        app,
        ["-g", "ATGCATGCATGCATGCATGCAT", "-b", "pINDUCER"],
    )
    without_loop = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "--loop",
            "",
        ],
    )
    assert with_loop.exit_code == 0 and without_loop.exit_code == 0
    assert "TAGTGAAGCCACAGATGTA" in with_loop.stdout
    assert "TAGTGAAGCCACAGATGTA" not in without_loop.stdout


def test_oligo_rejects_unknown_format():
    """An unrecognized -f value exits 1 with an 'Unsupported format' error."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "-f",
            "bogus",
        ],
    )
    assert result.exit_code == 1
    assert "Unsupported format" in result.stdout


def test_oligo_rejects_guide_and_batch_together(tmp_path: Path):
    """-g and --batch are mutually exclusive; passing both exits 1 with a
    'not both' error."""
    csv_path = tmp_path / "batch.csv"
    csv_path.write_text("name,antisense\nsh-1,ATGCATGCATGCATGCATGCAT\n")
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "--batch",
            str(csv_path),
            "-b",
            "pINDUCER",
        ],
    )
    assert result.exit_code == 1
    assert "not both" in result.stdout


@pytest.mark.parametrize("raw", ["abc", "1,1", "1,99"])
def test_oligo_rejects_invalid_mismatch_spec(raw):
    """A malformed -m value (non-numeric, or a position pair that's equal or
    out of range) exits 1 mentioning --mismatch."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "-m",
            raw,
        ],
    )
    assert result.exit_code == 1
    assert "--mismatch" in result.stdout


def test_keep_sites_adds_the_recognition_bases_to_the_oligos():
    """--keep-sites makes the insert carry each site's far half.

    pLKO cut with AgeI: the forward oligo gains AgeI's trailing T after the
    CCGG overhang, so ligation rebuilds ACCGGT instead of scarring it.
    """
    args = ["-g", "ATGCATGCATGCATGCATGCAT", "-b", "pLKO", "-f", "fasta"]
    plain = runner.invoke(app, args)
    kept = runner.invoke(app, [*args, "--keep-sites"])
    assert plain.exit_code == 0 and kept.exit_code == 0

    def forward(result):
        lines = [ln for ln in result.stdout.splitlines() if ln and not ln.startswith(">")]
        return lines[0]

    plain_fwd, kept_fwd = forward(plain), forward(kept)
    # AgeI's trailing T goes between the CCGG overhang and the insert body;
    # EcoRI's leading G goes on the far end, before the backbone's AATTC.
    assert plain_fwd.startswith("CCGG")
    assert kept_fwd == "CCGGT" + plain_fwd[len("CCGG") :] + "G"


def test_keep_sites_defaults_to_the_backbones_own_setting():
    """The flag is tri-state: unset means "whatever the backbone stores".

    pINDUCER ships with keep_sites=true, so a plain run already restores its
    XhoI site — and --no-keep-sites has to be able to turn that back off, or a
    stored True would be impossible to override.
    """
    args = ["-g", "TTGCACTGAAACTGGGTATCGG", "-b", "pINDUCER", "-f", "fasta"]

    def forward(result):
        """First record's sequence, rejoined.

        Rich wraps at the console width, and pINDUCER's oligos are ~119 nt —
        so the sequence arrives split across several lines, not one.
        """
        assert result.exit_code == 0, result.stdout
        lines = result.stdout.splitlines()
        start = next(i for i, ln in enumerate(lines) if ln.startswith(">")) + 1
        seq = []
        for line in lines[start:]:
            if line.startswith(">"):
                break
            seq.append(line.strip())
        return "".join(seq)

    default = forward(runner.invoke(app, args))
    explicit = forward(runner.invoke(app, [*args, "--keep-sites"]))
    disabled = forward(runner.invoke(app, [*args, "--no-keep-sites"]))

    # XhoI's G restored between the TCGA overhang and the flank, and EcoRI's
    # G on the far end — so disabling strips one base from each end.
    assert default.startswith("TCGAG") and default.endswith("G")
    assert default == explicit
    assert disabled == "TCGA" + default[len("TCGAG") : -1]


def test_keep_sites_warns_when_a_junction_has_no_site():
    """A directly stated overhang has no recognition sequence to rebuild, and
    must say so rather than let the flag quietly do nothing at that end."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "--site-5",
            "5'TCGA",
            "--keep-sites",
        ],
    )
    assert result.exit_code == 0
    assert "Nothing to restore at the 5' junction" in result.stderr


def test_degenerate_site_outside_the_overhang_is_accepted():
    """`N/TCGAG` still pins the overhang to TCGA, so it must not be rejected —
    only a degenerate code landing *on* the overhang is fatal."""
    ok = runner.invoke(
        app,
        ["-g", "ATGCATGCATGCATGCATGCAT", "-b", "pINDUCER", "--site-5", "N/TCGAG"],
    )
    assert ok.exit_code == 0, ok.stdout

    bad = runner.invoke(
        app,
        ["-g", "ATGCATGCATGCATGCATGCAT", "-b", "pINDUCER", "--site-5", "C/NNGTA"],
    )
    assert bad.exit_code == 1
    assert "overhang itself must be concrete" in bad.stdout


def test_oligo_rejects_invalid_site_override():
    """A --site-5 sequence with no recognizable cut position exits 1 with
    'No cut position'."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "--site-5",
            "GAATTC",
        ],
    )
    assert result.exit_code == 1
    assert "No cut position" in result.stdout


def test_oligo_json_format_content(tmp_path: Path):
    """`-f json -o <path>` writes a one-element list whose entry carries the
    guide sequence plus `forward`/`reverse` oligo fields."""
    out = tmp_path / "oligo"
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pINDUCER",
            "-f",
            "json",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    payload = json.loads(out.with_suffix(".json").read_text())
    assert len(payload) == 1
    assert payload[0]["guide"] == "ATGCATGCATGCATGCATGCAT"
    assert "forward" in payload[0] and "reverse" in payload[0]


def test_oligo_rna_does_not_corrupt_backbone_name_with_t():
    """Regression: --rna used to do a whole-text T->U replace, corrupting
    backbone names containing a literal T (e.g. pTRIPZ)."""
    result = runner.invoke(
        app,
        [
            "-g",
            "ATGCATGCATGCATGCATGCAT",
            "-b",
            "pTRIPZ",
            "--rna",
        ],
    )
    assert result.exit_code == 0
    assert "pTRIPZ" in result.stdout
    assert "pURIPZ" not in result.stdout


def test_oligo_batch_file_not_found():
    """A --batch path that doesn't exist exits 1 with a 'not found' error."""
    result = runner.invoke(
        app, ["--batch", "/nonexistent/path.csv", "-b", "pINDUCER"]
    )
    assert result.exit_code == 1
    assert "not found" in result.stdout


def test_oligo_batch_empty_csv(tmp_path: Path):
    """A --batch CSV with a header but no data rows exits 1 with 'No rows found'."""
    csv_path = tmp_path / "empty.csv"
    csv_path.write_text("name,antisense\n")
    result = runner.invoke(app, ["--batch", str(csv_path), "-b", "pINDUCER"])
    assert result.exit_code == 1
    assert "No rows found" in result.stdout


def test_oligo_batch_missing_guide_column(tmp_path: Path):
    """A --batch CSV lacking any recognizable guide/antisense column exits 1
    with 'Could not find guide column'."""
    csv_path = tmp_path / "bad_columns.csv"
    csv_path.write_text("name,sequence\nsh-1,ATGCATGCATGCATGCATGCAT\n")
    result = runner.invoke(app, ["--batch", str(csv_path), "-b", "pINDUCER"])
    assert result.exit_code == 1
    assert "Could not find guide column" in result.stdout


def test_oligo_batch_all_rows_invalid(tmp_path: Path):
    """A --batch CSV where every row has a blank guide exits 1 with
    'No valid guide sequences found'."""
    csv_path = tmp_path / "blank_guides.csv"
    csv_path.write_text("name,antisense\nsh-1,\nsh-2,\n")
    result = runner.invoke(app, ["--batch", str(csv_path), "-b", "pINDUCER"])
    assert result.exit_code == 1
    assert "No valid guide sequences found" in result.stdout


def test_oligo_batch_rejects_mismatch_missing_from_a_shorter_guide(tmp_path: Path):
    """A batch mismatch must be valid for every user-supplied guide, not only
    the first row."""
    csv_path = tmp_path / "mixed_lengths.csv"
    csv_path.write_text(
        "name,antisense\nlong,ATGCATGCATGCATGCATGCATGCAT\nshort,ATGCATGCATGCATGCATG\n"
    )
    result = runner.invoke(
        app,
        [
            "--batch",
            str(csv_path),
            "-b",
            "pINDUCER",
            "--mismatch",
            "25",
        ],
    )
    assert result.exit_code == 1
    assert "--mismatch" in result.stdout


def test_oligo_guide_rejects_fasta_input_conflict(tmp_path: Path):
    """INPUT (a FASTA to score) and --guide/--batch design oligos in two
    fundamentally different ways; combining them exits 1 rather than
    silently ignoring one."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "-g", "ATGCATGCATGCATGCATGCAT", "-b", "pINDUCER"],
    )
    assert result.exit_code == 1
    assert "--guide/--batch" in result.stdout


# ---------------------------------------------------------------------------
# Auto-triggered oligo design attached to scoring (--backbone + INPUT)
# ---------------------------------------------------------------------------


def test_auto_oligo_design_prints_panel_in_table_output():
    """--backbone alongside a FASTA scores as usual and also prints a
    separate 'Oligo Pairs' panel below the candidates table."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "-t", "2", "--multitarget", "1", "-b", "pINDUCER", "-f", "table"],
    )
    assert result.exit_code == 0
    assert "Oligo Pairs" in result.stdout
    assert "Forward" in result.stdout and "Reverse" in result.stdout


def test_auto_oligo_design_csv_writes_a_second_oligos_file(tmp_path: Path):
    """CSV can't hold two tables, so auto-designed oligos land in a second
    `*_oligos.csv` file next to the candidates CSV."""
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "2",
            "--multitarget",
            "1",
            "-b",
            "pINDUCER",
            "-f",
            "csv",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    candidates_csv = out.with_suffix(".csv")
    oligos_csv = out.parent / f"{out.name}_oligos.csv"
    assert candidates_csv.exists()
    assert oligos_csv.exists()
    text = oligos_csv.read_text()
    assert "forward" in text and "reverse" in text


def test_auto_oligo_design_json_embeds_oligos_key(tmp_path: Path):
    """JSON can hold two sections in one file: an 'oligos' key sits
    alongside 'candidates', no second file needed."""
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "2",
            "--multitarget",
            "1",
            "-b",
            "pINDUCER",
            "-f",
            "json",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    data = json.loads(out.with_suffix(".json").read_text())
    assert "candidates" in data and "oligos" in data
    assert len(data["oligos"]) == 2
    assert "forward" in data["oligos"][0] and "reverse" in data["oligos"][0]


def test_scoring_attached_oligo_design_uses_the_backbones_mismatch_position(
    tmp_path: Path,
):
    """Oligo design auto-triggered by --backbone (no -m) picks up that
    backbone's own `mismatch_positions`, and scoring's passenger preview
    agrees with it — same invariant as the explicit-flag case, just sourced
    from the backbone instead."""
    save_backbone(
        Backbone(
            name="Scoring Position Two Backbone",
            five_prime_enzyme="",
            five_prime_overhang=parse_overhang("G/AATTC"),
            five_prime_flank="AAA",
            loop="TAGTGAAGCCACAGATGTA",
            three_prime_flank="CCC",
            three_prime_enzyme="",
            three_prime_overhang=parse_overhang("C/TCGAG"),
            guide_on_5p=False,
            mismatch_positions=[2],
        )
    )
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "1",
            "--multitarget",
            "1",
            "-b",
            "Scoring Position Two Backbone",
            "-f",
            "json",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    data = json.loads(out.with_suffix(".json").read_text())

    guide = data["candidates"][0]["guide_sequence"]
    rc = str(Seq(guide).reverse_complement())
    scoring_passenger = data["candidates"][0]["passenger_sequence"]
    oligo_passenger = data["oligos"][0]["passenger"]

    for passenger in (scoring_passenger, oligo_passenger):
        assert passenger[0] == rc[0]
        assert passenger[1] != rc[1]


def test_auto_oligo_design_xlsx_adds_a_second_sheet(tmp_path: Path):
    """XLSX gains a second 'Oligo Pairs' sheet in the same workbook."""
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "2",
            "--multitarget",
            "1",
            "-b",
            "pINDUCER",
            "-f",
            "xlsx",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    import openpyxl

    wb = openpyxl.load_workbook(out.with_suffix(".xlsx"))
    assert "Oligo Pairs" in wb.sheetnames


def test_auto_oligo_design_fasta_appends_oligo_records():
    """FASTA can't hold two tables either, but it can hold more records —
    oligo fwd/rev records are appended after the guide records."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "-t", "2", "--multitarget", "1", "-b", "pINDUCER", "-f", "fasta"],
    )
    assert result.exit_code == 0
    assert "fwd_pINDUCER" in result.stdout
    assert "rev_pINDUCER" in result.stdout


def test_auto_oligo_design_manual_insert_flags_also_trigger_it():
    """A manual --loop/--adapter-5/--adapter-3 (no --backbone) triggers
    auto oligo design exactly like --backbone does."""
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "1",
            "--multitarget",
            "1",
            "--loop",
            "TAGTGAAGCCACAGATGTA",
            "--adapter-5",
            "",
            "--adapter-3",
            "",
            "-f",
            "table",
        ],
    )
    assert result.exit_code == 0
    assert "Oligo Pairs" in result.stdout


def test_top_n_default_is_ten_without_oligo_design():
    """Plain scoring (no --backbone/manual insert flags) keeps the
    pre-existing default of the top 10 candidates."""
    result = runner.invoke(
        app, [str(FIXTURE), "--multitarget", "1", "-f", "fasta"]
    )
    assert result.exit_code == 0
    n = sum(1 for line in result.stdout.splitlines() if line.startswith(">"))
    assert n == 10


def test_top_n_default_is_unlimited_with_oligo_design():
    """Once --backbone triggers oligo design, the default becomes every
    passing candidate rather than 10 — designing oligos for only a tenth
    of what passed filters would be a silent, surprising cap."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "--multitarget", "1", "-b", "pINDUCER", "-f", "fasta"],
    )
    assert result.exit_code == 0
    n_guides = sum(
        1
        for line in result.stdout.splitlines()
        if line.startswith(">") and "_fwd_" not in line and "_rev_" not in line
    )
    assert n_guides > 10


def test_top_n_explicit_still_caps_oligo_design(tmp_path: Path):
    """-t explicitly given always wins, even with oligo design triggered —
    checked via the oligo CSV's row count, which is unambiguous."""
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            str(FIXTURE),
            "-t",
            "1",
            "--multitarget",
            "1",
            "-b",
            "pINDUCER",
            "-f",
            "csv",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    oligos_csv = out.parent / f"{out.name}_oligos.csv"
    rows = [ln for ln in oligos_csv.read_text().splitlines() if ln]
    assert len(rows) == 2  # header + exactly one oligo pair


def test_scoring_mismatch_accepts_multiple_positions():
    """The merged --mismatch flag now accepts comma/space-separated multiple
    positions for scoring too, not just a single int."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "-t", "0", "-f", "table", "--mismatch", "1,3"],
    )
    assert result.exit_code == 0


def test_scoring_mismatch_zero_disables_mismatch():
    """--mismatch 0 is accepted and disables mismatch introduction, rather
    than being rejected as an out-of-range position."""
    result = runner.invoke(
        app,
        [str(FIXTURE), "-t", "0", "-f", "table", "--mismatch", "0"],
    )
    assert result.exit_code == 0


def test_oligo_mismatch_zero_yields_perfect_complement_passenger(tmp_path: Path):
    """-m 0 produces a passenger strand that's the guide's exact reverse
    complement, with no seed-region substitution."""
    guide = "ATGCATGCATGCATGCATGCAT"
    out = tmp_path / "oligo"
    result = runner.invoke(
        app,
        ["-g", guide, "-b", "pINDUCER", "-m", "0", "-f", "json", "-o", str(out)],
    )
    assert result.exit_code == 0
    payload = json.loads(out.with_suffix(".json").read_text())
    expected_passenger = str(Seq(guide).reverse_complement())
    assert payload[0]["passenger"] == expected_passenger


def _save_custom_backbone(name: str, mismatch_positions: list[int]) -> None:
    """A minimal blunt-ended backbone for tests exercising the
    `mismatch_positions` CSV column, without pulling in a real enzyme site."""
    save_backbone(
        Backbone(
            name=name,
            five_prime_enzyme="",
            five_prime_overhang=parse_overhang("G/AATTC"),
            five_prime_flank="AAA",
            loop="TAGTGAAGCCACAGATGTA",
            three_prime_flank="CCC",
            three_prime_enzyme="",
            three_prime_overhang=parse_overhang("C/TCGAG"),
            guide_on_5p=False,
            mismatch_positions=mismatch_positions,
        )
    )


def test_oligo_uses_the_backbones_own_mismatch_position_by_default(tmp_path: Path):
    """Without -m, a --backbone's stored `mismatch_positions` applies —
    not the flag's own "1" default."""
    guide = "ATGCATGCATGCATGCATGCAT"
    _save_custom_backbone("Position Two Backbone", [2])
    out = tmp_path / "oligo"
    result = runner.invoke(
        app,
        ["-g", guide, "-b", "Position Two Backbone", "-f", "json", "-o", str(out)],
    )
    assert result.exit_code == 0
    payload = json.loads(out.with_suffix(".json").read_text())

    rc = str(Seq(guide).reverse_complement())
    passenger = payload[0]["passenger"]
    # Position 1 (index 0) must be untouched; position 2 (index 1) mutated.
    assert passenger[0] == rc[0]
    assert passenger[1] != rc[1]


def test_explicit_mismatch_overrides_the_backbones_default(tmp_path: Path):
    """An explicit -m always wins over the backbone's stored default."""
    guide = "ATGCATGCATGCATGCATGCAT"
    _save_custom_backbone("Position Two Backbone Override", [2])
    out = tmp_path / "oligo"
    result = runner.invoke(
        app,
        [
            "-g",
            guide,
            "-b",
            "Position Two Backbone Override",
            "-m",
            "1",
            "-f",
            "json",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0
    payload = json.loads(out.with_suffix(".json").read_text())

    rc = str(Seq(guide).reverse_complement())
    passenger = payload[0]["passenger"]
    # -m 1 wins: position 1 mutated, position 2 (the backbone's own default) untouched.
    assert passenger[0] != rc[0]
    assert passenger[1] == rc[1]


def test_backbone_mismatch_position_out_of_range_for_guide_is_reported(
    tmp_path: Path,
):
    """A backbone whose stored mismatch position doesn't fit the guide being
    designed exits 1 with a clear diagnostic, rather than silently doing
    nothing."""
    _save_custom_backbone("Too Long Mismatch Backbone", [99])
    result = runner.invoke(
        app,
        ["-g", "ATGCATGCATGCATGCATGCAT", "-b", "Too Long Mismatch Backbone"],
    )
    assert result.exit_code == 1
    assert "Too Long Mismatch Backbone" in result.stdout
    assert "--mismatch" in result.stdout


# ---------------------------------------------------------------------------
# Bare invocation
# ---------------------------------------------------------------------------


def test_bare_shrna_invokes_picker(monkeypatch):
    """Bare ``knockknock`` should call into the TUI launch path (with no candidates)
    rather than printing CLI help. We intercept the call so no actual terminal
    UI is spawned during the test."""
    calls = []

    def fake_launch(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr("knockknock.cli.tui_launch", fake_launch)
    result = runner.invoke(app, [])
    assert result.exit_code == 0
    assert len(calls) == 1
    # Phase 3.1 contract: called with no candidates (None) → picker mode.
    args, kwargs = calls[0]
    assert (args and args[0] is None) or kwargs.get("candidates") is None


# ---------------------------------------------------------------------------
# Default-command routing: `run` is optional
# ---------------------------------------------------------------------------


def test_bare_positional_input_routes_to_run_without_the_word_run():
    """`knockknock targets.fasta` (no `run`) scores exactly like
    `knockknock run targets.fasta` — this is the whole point of the merge."""
    with_run = runner.invoke(app, ["run", str(FIXTURE), "-t", "5", "-f", "table"])
    without_run = runner.invoke(app, [str(FIXTURE), "-t", "5", "-f", "table"])
    assert with_run.exit_code == 0 and without_run.exit_code == 0
    assert _passed(with_run.stdout) == _passed(without_run.stdout)


def test_bare_option_first_still_routes_to_run():
    """An option before the positional (no `run`) must still route
    correctly — Click's own subcommand resolution would otherwise choke on
    a leading option it doesn't recognize at the group level."""
    result = runner.invoke(app, ["-f", "table", "-t", "5", str(FIXTURE)])
    assert result.exit_code == 0
    assert "Run Summary" in result.stdout


def test_bare_dash_stdin_routes_to_run():
    """A bare `-` (no `run`) still means "read stdin", not "unknown command"."""
    result = runner.invoke(
        app, ["-", "-t", "0"], input=">x\nATGCATGCATGCATGCATGCATGC\n"
    )
    assert result.exit_code == 0


def test_bare_oligo_flags_route_to_run():
    """`--guide`/`--backbone` alone (no `run`, no INPUT) still reach the
    standalone oligo design path."""
    result = runner.invoke(
        app, ["--guide", "ATGCATGCATGCATGCATGCAT", "-b", "pINDUCER"]
    )
    assert result.exit_code == 0
    assert "_fwd_" in result.stdout and "_rev_" in result.stdout


def test_bare_word_info_is_treated_as_a_filename_not_a_command():
    """`info` is a `--info` flag now, not a subcommand — the bare word
    `info` (no such file) is treated as INPUT and reported as missing,
    confirming there's no special-cased command name left to collide with."""
    result = runner.invoke(app, ["info"])
    assert result.exit_code == 1
    assert "not found" in result.stdout


def test_run_is_hidden_from_top_level_help():
    """`run` never needs to be typed, so it shouldn't clutter `--help`'s
    command list — the full option reference is still reachable via
    `knockknock run --help` for anyone who wants it."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "--info" in result.stdout  # the merged help text mentions it
    assert "Commands" not in result.stdout
