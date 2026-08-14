"""End-to-end tests that exec the CLI as a real, separate OS process.

Every other CLI test in this package (`test_cli.py`) drives `knockknock.cli.app`
in-process via Click's `CliRunner`. That's fast and covers argument parsing/
business logic thoroughly, but `CliRunner` fakes stdin as an in-memory buffer
with no real file descriptor. `cli.py`'s bare-invocation TTY/pipe detection
(`_entry`) has an explicit carve-out for exactly that case: `sys.stdin.fileno()`
raises under `CliRunner`, and the code treats that as "must be a test
environment" and takes the TUI-launch branch rather than the real pipe-detection
branch. That means the *actual* code path a real `cat targets.fa | kk` hits —
a genuine non-tty stdin with a real fd, auto-scoring with no subcommand and no
`-` — is structurally untestable through `CliRunner`. These tests exec a real
subprocess instead, so stdin/stdout/stderr and exit codes are the same OS-level
things an actual user's shell would see.

Each subprocess pays Python's own startup cost plus this package's imports
(pandas/openpyxl/textual, ...), so keep this file to a handful of targeted
smoke tests — the exhaustive option/behavior matrix belongs in `test_cli.py`,
not here.
"""

import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "data" / "sample.fa"


def _run(args: list[str], input: str | None = None) -> subprocess.CompletedProcess:
    """Run the CLI in a fresh interpreter, the same way the installed `kk`/
    `knockknock` console scripts do (`from knockknock.cli import app;
    sys.exit(app())`) — without depending on those scripts actually being on
    PATH, which isn't guaranteed in every environment `pytest` runs in."""
    return subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from knockknock.cli import app; sys.exit(app(prog_name='kk'))",
            *args,
        ],
        input=input,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_help_runs_as_a_real_process():
    result = _run(["--help"])
    assert result.returncode == 0
    assert "Usage: kk" in result.stdout
    assert "--backbone" in result.stdout
    assert "Traceback" not in result.stderr


def test_scoring_smoke_as_a_real_process():
    result = _run([str(FIXTURE), "-t", "5", "-f", "table"])
    assert result.returncode == 0
    assert "Passed all filters" in result.stdout


def test_real_piped_stdin_is_auto_scored_with_no_subcommand():
    """`cat targets.fa | kk` (bare — no `run`, no `-`) on a real OS pipe:
    non-tty stdin with a genuine file descriptor should hit `_entry`'s
    auto-score branch, not the TUI-launch branch. This is the one behavior
    `CliRunner` cannot exercise at all (see module docstring).

    A single raw (unheadered) sequence, not the multi-record fixture:
    `_entry`'s auto-score defaults `multitarget` to "present in every input
    mRNA", so with `FIXTURE`'s 4 unrelated records almost nothing passes and
    stdout would be empty regardless of which branch ran — not a useful
    assertion for "did this hit auto-score at all".
    """
    result = _run([], input="ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
    assert result.returncode == 0
    assert result.stdout.startswith(">")


def test_real_piped_stdin_dash_reads_as_fasta():
    result = _run(["-", "-t", "0"], input=FIXTURE.read_text())
    assert result.returncode == 0


def test_csv_export_writes_a_real_file(tmp_path: Path):
    out = tmp_path / "out"
    result = _run([str(FIXTURE), "-t", "3", "-f", "csv", "-o", str(out)])
    assert result.returncode == 0
    assert out.with_suffix(".csv").exists()


def test_missing_input_file_exits_nonzero_with_a_message():
    result = _run([str(FIXTURE.parent / "does_not_exist.fasta")])
    assert result.returncode != 0
    assert "Traceback" not in result.stderr
