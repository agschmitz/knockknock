"""Typer CLI — a single merged command plus `info`."""

import io
import os
import signal
import sys
import tempfile
from importlib.resources import files
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table
from typer import rich_utils
from typer.core import HAS_RICH, TyperGroup

# Python installs a SIGPIPE handler that turns a closed downstream pipe (e.g.
# `kk targets.fasta | head`) into a BrokenPipeError instead of just killing us
# — restore the default disposition so piping into a short-lived reader exits
# cleanly, the way `grep`/`git` etc. do. No-op on Windows, which has no SIGPIPE.
if hasattr(signal, "SIGPIPE"):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

from knockknock import __version__
from knockknock.backbone import (
    Backbone,
    backbones_data_path,
    bundled_backbones_path,
    get_backbone,
    list_backbone_names,
    load_backbones,
    restore_bundled_backbones,
)
from knockknock.cutsite import (
    BLUNT,
    LEFT,
    RIGHT,
    CutSiteError,
    Overhang,
    invalid_dna_chars,
    parse_overhang,
    validate_dna_sequence,
)
from knockknock.engine import (
    DEFAULT_GUIDE_LENGTH,
    MAX_GUIDE_LENGTH,
    MIN_GUIDE_LENGTH,
    QUALITY_PRESETS,
    AnalysisParameterError,
    AnalysisParams,
    make_shrnas,
    process_shrnas,
    validate_analysis_params,
    validate_top_n,
)
from knockknock.export import (
    no_candidates_message,
    resolve_export_path,
    to_csv_text,
    to_fasta,
    to_json_text,
    to_terminal,
    to_xlsx,
)
from knockknock.oligo import (
    BatchInputError,
    OligoConfig,
    OligoPair,
    assemble_batch,
    assemble_oligos,
    check_ends,
    check_insert,
    check_sites,
    export_oligos_csv_text,
    export_oligos_fasta_text,
    export_oligos_json_text,
    parse_mismatch_positions,
    read_batch_csv,
)
from knockknock.reference import load_coeffs, load_deltaG, load_dG_bins
from knockknock.seqinput import (
    LoadedSequences,
    SequenceInputError,
    load_sequences,
    read_sequence_file,
)
from knockknock.tui import launch as tui_launch

console = Console()
err_console = Console(stderr=True)


# Long-form --help body. Kept as a plain-text resource file, not an inlined
# Python string, so it can be read (and edited) without touching this module —
# the CLI's own imports and error handling are the only thing that belongs
# here. `\b` markers inside it are Click's "don't rewrap this block" escape,
# preserved verbatim by reading the file as-is.
_APP_HELP = files("knockknock").joinpath("cli_help.txt").read_text(encoding="utf-8")


class _DefaultCommandGroup(TyperGroup):
    """Click group that falls back to the ``run`` command when the first
    token isn't a registered subcommand name.

    This is what lets ``kk targets.fasta ...`` and ``kk --backbone ...
    targets.fasta`` work without spelling out ``run`` — Click otherwise
    requires an exact registered command name as the first token. ``run``
    is the only command ever registered (``info`` is a ``--info`` flag, not
    a subcommand — see `run`'s docstring), and it's kept `hidden=True`
    since it never needs to be typed; explicit ``kk run targets.fasta``
    still works unchanged. Pairs with ``ignore_unknown_options=True`` on the
    app's context so a leading option (e.g. ``kk --backbone ...``) doesn't
    get rejected by the group's own (option-free) argument parser before
    reaching this fallback.
    """

    default_command = "run"

    def resolve_command(self, ctx, args):
        if args and args[0] not in self.commands:
            args = [self.default_command, *args]
        return super().resolve_command(ctx, args)

    def format_help(self, ctx, formatter):
        """Render ``kk --help`` as if it were ``kk run --help``.

        Click/Typer never fold a subcommand's Arguments/Options panels into
        its group's own help — a group's help only ever reflects its *own*
        params, which here is just ``--help`` since every real option lives
        on ``run``. Since ``run`` is deliberately hidden (nobody should have
        to type it — see the class docstring) but is where every scoring/
        oligo-design flag actually lives, the naive result is a `--help`
        that documents nothing. ``ctx`` here is already the group's own root
        context (``ctx.command_path`` is just ``kk``, no parent) — reusing
        it with ``run``'s command swapped in is what makes typer's own
        `rich_format_help` render ``run``'s Arguments/Options under the
        group's program name, with no extra context object needed.
        """
        if not (HAS_RICH and self.rich_markup_mode is not None):
            return super().format_help(ctx, formatter)
        run_cmd = self.commands[self.default_command]
        rich_utils.rich_format_help(
            obj=run_cmd, ctx=ctx, markup_mode=self.rich_markup_mode
        )


app = typer.Typer(
    name="knockknock",
    cls=_DefaultCommandGroup,
    context_settings={
        "ignore_unknown_options": True,
        # "-h" as a short alias for "--help". Click inherits
        # `help_option_names` down to child contexts (here, `run`'s) when
        # the child doesn't set its own — see Context.__init__ — so this one
        # place covers both `kk -h` and `kk run -h`.
        "help_option_names": ["-h", "--help"],
    },
    help=_APP_HELP,
    add_completion=False,
)


@app.callback(invoke_without_command=True)
def _entry(ctx: typer.Context) -> None:
    """When invoked with no subcommand, auto-score piped input or open the TUI."""
    if ctx.invoked_subcommand is None:
        # First detect TTY: if `fileno()` works and `isatty()` is True, we're
        # in an interactive terminal → TUI.  Otherwise (pipe, file redirect,
        # or a test environment where `fileno()` raises), check whether
        # `fileno()` actually gave us a real fd — if it raised, it's pytest
        # → TUI; otherwise it's real data → auto-score.
        try:
            is_tty = sys.stdin.isatty()
        except (OSError, ValueError, io.UnsupportedOperation):
            # No real fd (e.g. pytest) — TUI path.
            is_tty = True
        else:
            if not is_tty:
                # Not a TTY — but is it a real pipe/redirect or a test env
                # that mimics one?  Test env `fileno()` raises; real pipe
                # succeeds.
                try:
                    sys.stdin.fileno()
                except (OSError, ValueError, io.UnsupportedOperation):
                    # No fd → test env → TUI
                    is_tty = True
                # else — fd exists and not a TTY → `kk` or `cat | kk` → auto-score
        if is_tty:
            tui_launch(None)
            raise typer.Exit()
        # Auto-score: real non-TTY stdin (pipe / file redirect).
        try:
            text = _read_stdin_text()
        except typer.Exit:
            return
        try:
            loaded = load_sequences(text)
        except SequenceInputError as exc:
            console.print(f"[red]Error:[/] {exc}")
            raise typer.Exit(1) from exc
        mrna_tuples, input_warnings = loaded
        n_mrnas = len(mrna_tuples)
        # stderr so warnings don't contaminate piped FASTA output.
        for message in input_warnings:
            err_console.print(f"[yellow]Warning:[/] {message}")
        # Build minimal params with "strict" defaults and no filters.
        params = AnalysisParams(
            length=DEFAULT_GUIDE_LENGTH,
            gc_min=0.0,
            gc_max=100.0,
            min_efficacy_score=None,
            remove_motifs=False,
            multitarget=max(n_mrnas, 1),
            rna_output=False,
            stem_mismatch=[1],
        )
        candidates_by_mrna = make_shrnas(mrna_tuples, params)
        results, stats = process_shrnas(
            candidates_by_mrna, mrna_tuples, params
        )
        exported = results  # unlimited: output everything on a pipe
        if not results:
            err_console.print(
                f"[yellow]{no_candidates_message(stats, params)}[/yellow]"
            )
        sys.stdout.write(to_fasta(exported, params))
        sys.stdout.flush()
        return


def _normalise_guide(guide: str, *, label: str, exact_length: bool = False) -> str:
    """Validate a user-supplied guide and return its DNA representation."""
    cleaned = "".join(guide.split()).upper().replace("U", "T")
    if not cleaned:
        raise ValueError(f"{label} is empty.")
    degenerate, other = invalid_dna_chars(cleaned)
    invalid = "".join(sorted(set(degenerate) | set(other)))
    if invalid:
        raise ValueError(
            f"{label} contains invalid base(s) {invalid!r}; expected A, C, G, T, or U."
        )
    if exact_length and not MIN_GUIDE_LENGTH <= len(cleaned) <= MAX_GUIDE_LENGTH:
        raise ValueError(
            f"{label} must be {MIN_GUIDE_LENGTH}–{MAX_GUIDE_LENGTH} nt for this "
            f"scoring model; got {len(cleaned)}."
        )
    if not exact_length and not 19 <= len(cleaned) <= 30:
        raise ValueError(
            f"{label} must be 19–30 nt for oligo design; got {len(cleaned)}."
        )
    return cleaned


def _parse_mismatch_positions(raw: str, guide_length: int) -> list[int]:
    """CLI wrapper over `oligo.parse_mismatch_positions` that prefixes the
    generic validation message with the flag name that's actually wrong."""
    try:
        return parse_mismatch_positions(raw, guide_length)
    except ValueError as exc:
        raise ValueError(f"--mismatch {exc}") from exc


def _mismatch_param_is_explicit(ctx: typer.Context) -> bool:
    """Whether the user actually typed ``-m``/``--mismatch``, as opposed to
    it sitting at its ``"1"`` default.

    Needed to decide whether a ``--backbone``'s own `Backbone
    .mismatch_positions` column should apply: an explicit ``-m`` always
    wins, but the flag's own default must not silently shadow what the
    chosen backbone declares.
    """
    try:
        source = ctx.get_parameter_source("mismatch")
    except (AttributeError, KeyError):
        # Older/newer click without the API: treat as explicit, so this
        # falls back to the pre-existing (flag-only) behaviour rather than
        # a backbone default it can't actually confirm was unrequested.
        return True
    return not (source is None or source.name == "DEFAULT")


def _resolve_mismatch_positions(
    ctx: typer.Context,
    mismatch: str,
    guide_length: int,
    backbone_name: str | None,
) -> list[int]:
    """Resolve ``--mismatch`` to a concrete, guide-length-checked position
    list, letting a named backbone's own default apply when ``-m`` was left
    unset.

    An explicit ``-m``/``--mismatch`` always wins — parsed and range-checked
    exactly as before. Left at its default, a named ``--backbone``'s own
    `Backbone.mismatch_positions` is used instead, so a vector that's
    conventionally opened at a different seed position doesn't need ``-m``
    re-typed on every run; that list is still range-checked here (against
    *this* guide length, not the backbone-agnostic one `backbone.py` itself
    can check) so a stale/miscopied database row is caught immediately
    rather than silently doing nothing (`compute_mismatched_passenger`
    otherwise ignores out-of-range positions rather than raising). A
    backbone name that doesn't resolve is left for `_resolve_oligo_config`
    to report; this just falls through to the flag's own default.
    """
    if not _mismatch_param_is_explicit(ctx) and backbone_name:
        bb = get_backbone(backbone_name)
        if bb is not None:
            positions = list(bb.mismatch_positions)
            bad = [p for p in positions if not 1 <= p <= guide_length]
            if bad:
                raise ValueError(
                    f"backbone {bb.name!r} has mismatch_positions {bad} that "
                    f"don't fit a {guide_length}-nt guide; override with "
                    "--mismatch."
                )
            return positions
    return _parse_mismatch_positions(mismatch, guide_length)


def _read_stdin_text() -> str:
    """Read raw sequence text from stdin, prompting at an interactive terminal.

    Parsing and validation are not done here — they belong to
    `seqinput.load_sequences`, which every input route shares.
    """
    if sys.stdin.isatty():
        err_console.print(
            "[dim]Reading from stdin — paste a raw sequence or FASTA, then "
            "press Ctrl+D (Ctrl+Z then Enter on Windows) to finish.[/dim]"
        )
    raw_text = sys.stdin.read()
    if not raw_text.strip():
        console.print("[red]Error:[/] No input provided on stdin.")
        raise typer.Exit(1)
    return raw_text


def _load_input(fasta_input: str) -> LoadedSequences:
    """Turn the INPUT argument into validated, engine-ready mRNAs.

    ``-`` reads stdin; anything else is a file path. Both go through the same
    `seqinput.load_sequences`, so a FASTA behaves identically whether it was
    piped, pasted, or named on the command line.

    Non-fatal notes come back with the sequences rather than being printed
    here: `--tui` hands them to the TUI as toasts instead, since anything
    written to the terminal is wiped the moment the app starts.
    """
    if fasta_input == "-":
        text = _read_stdin_text()
    else:
        input_path = Path(fasta_input)
        if not input_path.is_file():
            console.print(f"[red]Error:[/] File '{fasta_input}' not found.")
            raise typer.Exit(1)
        try:
            text = read_sequence_file(input_path)
        except SequenceInputError as exc:
            console.print(f"[red]Error:[/] {exc}")
            raise typer.Exit(1) from exc

    try:
        return load_sequences(text)
    except SequenceInputError as exc:
        console.print(f"[red]Error:[/] {exc}")
        raise typer.Exit(1) from exc


def _resolve_output_path(user_output: str | None, fasta_input: str, ext: str) -> Path:
    return resolve_export_path(ext, user_output=user_output, reference_path=fasta_input)


def _oligo_csv_path(user_output: str | None, fasta_input: str) -> Path:
    """Where the auto-designed oligo CSV goes, alongside the candidates CSV.

    CSV can't hold two distinct tables cleanly in one file, so it's the one
    format that gets a second *file* rather than a second section — named
    from the same base as the candidates export (explicit ``-o`` or
    auto-named from the input) plus an ``_oligos`` tag.
    """
    if user_output and user_output != "-":
        return Path(f"{user_output}_oligos").expanduser().with_suffix(".csv").resolve()
    return resolve_export_path(".csv", reference_path=fasta_input, label="oligos")


def _dedupe_export_path(out_path: Path, seen: set[Path]) -> Path:
    """Disambiguate ``out_path`` against paths already written this run.

    Auto-named export paths are timestamped to the second, so a batch of
    inputs sharing a stem (or the same file listed twice) can resolve to the
    same path within one run and silently overwrite each other. Only
    renames on an actual collision — the common case (distinct stems, or
    any gap of a second or more) is untouched.
    """
    if out_path not in seen:
        seen.add(out_path)
        return out_path

    n = 2
    candidate = out_path.with_stem(f"{out_path.stem}_{n}")
    while candidate in seen:
        n += 1
        candidate = out_path.with_stem(f"{out_path.stem}_{n}")
    err_console.print(
        f"[yellow]Warning:[/] Output name collided with an earlier file in "
        f"this run; writing to {candidate} instead."
    )
    seen.add(candidate)
    return candidate


def _resolve_params(
    ctx: typer.Context,
    *,
    n_mrnas: int,
    length: int,
    quality: str,
    gc_min: float,
    gc_max: float,
    min_efficacy_score: float | None,
    sort_dg: bool,
    remove_motifs: bool,
    top_n: int | None,
    multitarget: int | None,
    mismatch: list[int],
    rna: bool,
    seed: int | None,
) -> AnalysisParams:
    """Apply the quality preset, range-check everything, and build AnalysisParams.

    Split out of `run` so the preset-versus-explicit-flag precedence rules and
    the numeric bounds are one readable block rather than eighty lines wedged
    between input parsing and output dispatch.
    """
    if quality not in QUALITY_PRESETS:
        console.print(
            f"[red]Error:[/] Unknown --quality '{quality}'. "
            f"Choose one of: {', '.join(QUALITY_PRESETS)}."
        )
        raise typer.Exit(1)

    # Apply the preset, but only to filter flags the user did NOT explicitly
    # set — click's parameter-source machinery is what tells us which values
    # actually came from the command line, as opposed to being the declared
    # default that happens to equal what the user wanted.
    preset = QUALITY_PRESETS[quality]

    def _from_default(name: str) -> bool:
        try:
            src = ctx.get_parameter_source(name)
        except (AttributeError, KeyError):
            # Older/newer click without the API, or a name we don't own.
            # Treat as "not set by the user" so the preset still applies.
            return True
        return src is None or src.name == "DEFAULT"

    if _from_default("gc_min"):
        gc_min = preset["gc_min"]
    if _from_default("gc_max"):
        gc_max = preset["gc_max"]
    if _from_default("remove_motifs"):
        remove_motifs = preset["remove_motifs"]
    if _from_default("min_efficacy_score") and preset["min_efficacy_score"] is not None:
        min_efficacy_score = preset["min_efficacy_score"]

    if multitarget is None:
        # Not passed on the command line: default to requiring the guide in
        # every input mRNA, same as the TUI's "target all inputs" default.
        multitarget = n_mrnas
    elif multitarget > n_mrnas:
        multitarget = n_mrnas
        console.print(
            f"[yellow]Warning:[/] multitarget reduced to {multitarget} "
            f"(number of input sequences)."
        )

    try:
        if top_n is not None:
            validate_top_n(top_n)
        return validate_analysis_params(
            length=length,
            gc_min=gc_min,
            gc_max=gc_max,
            min_efficacy_score=min_efficacy_score,
            sort_by_dg=sort_dg,
            remove_motifs=remove_motifs,
            multitarget=multitarget,
            rna_output=rna,
            mismatch=mismatch,
            seed=seed,
        )
    except AnalysisParameterError as exc:
        # The shared validator is frontend-neutral; preserve established CLI
        # diagnostics by adding the corresponding option name here.
        message = str(exc)
        replacements = {
            "length": "--length",
            "top": "--top",
            "multitarget": "--multitarget",
            "mismatch": "--mismatch",
        }
        for plain, option in replacements.items():
            if message.startswith(plain):
                message = option + message[len(plain) :]
                break
        if message.startswith("GC limits must satisfy"):
            message = "GC limits must satisfy 0 ≤ --gc-min < --gc-max ≤ 100."
        console.print(f"[red]Error:[/] {message}")
        raise typer.Exit(1) from exc


def _print_backbones_table() -> None:
    table = Table(title="Available Backbones")
    table.add_column("Name", style="cyan")
    table.add_column("5' Enzyme")
    table.add_column("5' End")
    table.add_column("5' Flank")
    table.add_column("Loop")
    table.add_column("3' Flank")
    table.add_column("3' Enzyme")
    table.add_column("3' End")
    table.add_column("Guide on 5'")
    table.add_column("Mismatch Pos")

    for b in load_backbones():
        table.add_row(
            b.name,
            b.five_prime_enzyme,
            b.five_prime_overhang.describe(),
            b.five_prime_flank,
            b.loop,
            b.three_prime_flank,
            b.three_prime_enzyme,
            b.three_prime_overhang.describe(),
            str(b.guide_on_5p),
            ",".join(str(m) for m in b.mismatch_positions) or "none",
        )
    # Sequence columns (enzyme sites, flanks, loop) are short but numerous;
    # the default 80-col terminal width truncates them into illegible
    # ellipses, so size the print to the table's unconstrained natural
    # width instead of the terminal width. `console.print(..., width=...)`
    # doesn't reliably widen a Table's own layout — setting `console.width`
    # directly (and restoring it after) does.
    unconstrained = console.options.update(max_width=10_000)
    natural_width = table.__rich_measure__(console, unconstrained).maximum
    saved_width = console.width
    console.width = max(console.width, natural_width)
    try:
        console.print(table)
    finally:
        console.width = saved_width


def _build_oligo_config(
    backbone: Backbone | None,
    *,
    loop: str | None,
    adapter_5: str | None,
    adapter_3: str | None,
    site_5: str | None,
    site_3: str | None,
    guide_on_5p: bool,
    keep_sites: bool | None = None,
) -> OligoConfig:
    """Start from ``backbone`` (or a bare blunt config) and apply the overrides.

    Overrides are tested with ``is not None``, not truthiness: an explicitly
    empty value (``--loop ''``) is a deliberate "this segment is absent" and
    must override a backbone's sequence, rather than being ignored as falsy.
    """
    if backbone is not None:
        # `mismatches` deliberately left unset: `_resolve_oligo_config`
        # always overwrites `config.mismatch_positions` with the caller's
        # resolved value right after this returns, so this only ever needs
        # to be *some* valid placeholder — `from_backbone`'s own default
        # (the backbone's stored `mismatch_positions`) is as good as any.
        config = OligoConfig.from_backbone(backbone, keep_sites=keep_sites)
    else:
        config = OligoConfig(
            five_prime_overhang=Overhang(kind=BLUNT, junction=LEFT),
            five_prime_flank="",
            loop="",
            three_prime_flank="",
            three_prime_overhang=Overhang(kind=BLUNT, junction=RIGHT),
            guide_on_5p=guide_on_5p,
            mismatch_positions=[1],
        )

    try:
        if loop is not None:
            config.loop = loop.strip().upper()
            validate_dna_sequence(config.loop, "Loop")
        if adapter_5 is not None:
            config.five_prime_flank = adapter_5.strip().upper()
            validate_dna_sequence(config.five_prime_flank, "5' flank")
        if adapter_3 is not None:
            config.three_prime_flank = adapter_3.strip().upper()
            validate_dna_sequence(config.three_prime_flank, "3' flank")
        if site_5 is not None:
            config.five_prime_overhang = parse_overhang(site_5, junction=LEFT)
        if site_3 is not None:
            config.three_prime_overhang = parse_overhang(site_3, junction=RIGHT)
    except CutSiteError as exc:
        console.print(f"[red]Error:[/] {exc}")
        raise typer.Exit(1) from exc

    if guide_on_5p:
        config.guide_on_5p = True
    if keep_sites is not None:
        config.keep_sites = keep_sites
    return config


def _resolve_oligo_config(
    backbone_name: str | None,
    *,
    loop: str | None,
    adapter_5: str | None,
    adapter_3: str | None,
    site_5: str | None,
    site_3: str | None,
    guide_on_5p: bool,
    mismatch_positions: list[int],
    keep_sites: bool | None = None,
) -> OligoConfig:
    """Resolve --backbone (or a full manual loop+adapter override) into an
    OligoConfig, or exit 1 with a diagnostic.

    Shared by both oligo design paths — the standalone ``--guide``/``--batch``
    mode and the auto-trigger attached to scoring — so a missing/unknown
    backbone is reported identically either way.
    """
    bb = None
    if backbone_name:
        bb = get_backbone(backbone_name)
        if bb is None:
            names = list_backbone_names()
            console.print(
                f"[red]Error:[/] Backbone '{backbone_name}' not found. "
                f"Available: {', '.join(names)}"
            )
            raise typer.Exit(1)

    # Without a named backbone there is no design strategy to fall back on:
    # the loop and both flanks must be stated explicitly. Empty strings are a
    # valid answer ("--loop ''" really does mean no loop) — what is rejected is
    # leaving them unsaid, which previously produced a loop-less, flank-less
    # guide-fused-to-passenger oligo that looked orderable. Defaulting to the
    # first row of the backbone database would be just as wrong: silently
    # picking someone else's vector is not a decision this tool should make.
    if bb is None:
        missing = [
            flag
            for flag, value in (
                ("--loop", loop),
                ("--adapter-5", adapter_5),
                ("--adapter-3", adapter_3),
            )
            if value is None
        ]
        if missing:
            names = list_backbone_names()
            console.print(
                "[red]Error:[/] No oligo design strategy given. Either name a "
                "backbone with --backbone, or define the insert yourself with "
                f"{', '.join(missing)} (missing).\n"
                f"  Available backbones: {', '.join(names) if names else '(none)'}\n"
                "  List them in full:    knockknock --list-backbones\n"
                "  Pass an empty value if a segment really is absent, "
                "e.g. --loop ''"
            )
            raise typer.Exit(1)

    config = _build_oligo_config(
        bb,
        loop=loop,
        adapter_5=adapter_5,
        adapter_3=adapter_3,
        site_5=site_5,
        site_3=site_3,
        guide_on_5p=guide_on_5p,
        keep_sites=keep_sites,
    )
    config.mismatch_positions = mismatch_positions
    return config


def _warn_about_design(config: OligoConfig, pairs: list[OligoPair]) -> None:
    """Emit design warnings to stderr, de-duplicated across the batch.

    stderr, not stdout, so warnings never contaminate piped sequence output.
    """
    for message in check_ends(config) + check_sites(config):
        err_console.print(f"[yellow]Warning:[/] {message}")

    seen: set[str] = set()
    for pair in pairs:
        arm1, arm2 = config.arms(pair.guide, pair.passenger)
        for message in check_insert(config, arm1, arm2):
            if message in seen:
                continue
            seen.add(message)
            err_console.print(f"[yellow]Warning ({pair.name}):[/] {message}")


def _emit_oligos(
    pairs: list[OligoPair],
    config: OligoConfig,
    *,
    output_format: str,
    output: str | None,
    rna: bool,
) -> None:
    """Render the assembled pairs and either write them or print them."""
    if output_format == "json":
        text = export_oligos_json_text(pairs, config.backbone_name, rna=rna)
        suffix = ".json"
    else:
        text = export_oligos_fasta_text(pairs, config.backbone_name, rna=rna)
        suffix = ".fasta"

    if not output:
        console.print(text)
        return

    out_path = Path(output).with_suffix(suffix)
    out_path.write_text(text, encoding="utf-8")
    console.print(f"[green]Exported {len(pairs)} pairs to {out_path}[/green]")


@app.command(name="run", hidden=True, help=_APP_HELP)
def run(
    ctx: typer.Context,
    fasta_inputs: list[str] | None = typer.Argument(
        None,
        metavar="INPUT...",
        help="One or more FASTA files with mRNA sequences to score. Use '-', "
        "or omit entirely, to read from stdin. Must be omitted when using "
        "--guide/--batch, which design oligos directly instead of scoring "
        "anything.",
    ),
    output: str | None = typer.Option(
        None, "-o", "--output", help="Output file path (without extension)"
    ),
    output_format: str | None = typer.Option(
        None,
        "-f",
        "--format",
        help="Scoring mode: table, fasta, csv, xlsx, json (default 'table' at "
        "a terminal, 'fasta' when piped). --guide/--batch mode: fasta, json "
        "(default 'fasta').",
    ),
    quality: str = typer.Option(
        "strict",
        "--quality",
        help="Filter preset: lax | standard | strict. Individual flags override.",
    ),
    length: int = typer.Option(
        DEFAULT_GUIDE_LENGTH,
        "-l",
        "--length",
        help=f"Guide strand length in nt ({MIN_GUIDE_LENGTH}–{MAX_GUIDE_LENGTH}). "
        "Efficacy scores are sums over positions, so they are only "
        "comparable between candidates of the same length.",
    ),
    gc_min: float = typer.Option(0.0, "--gc-min", help="Minimum GC percentage (0–100)"),
    gc_max: float = typer.Option(100.0, "--gc-max", help="Maximum GC percentage (0–100)"),
    min_efficacy_score: float | None = typer.Option(
        None,
        "-s",
        "--min-efficacy-score",
        help="Minimum efficacy score (no minimum if unset)",
    ),
    sort_dg: bool = typer.Option(
        False,
        "--sort-dg",
        help="Sort by ΔG score instead of efficacy score",
    ),
    remove_motifs: bool = typer.Option(
        False,
        "--remove-motifs",
        help="Filter out homopolymeric and low-complexity motifs",
    ),
    top_n: int | None = typer.Option(
        None,
        "-t",
        "--top",
        help="Limit to the top N candidates — and, when oligos are designed, "
        "the N that get them. Default: 10 for plain scoring; unlimited (every "
        "passing candidate) the moment oligo design is triggered.",
    ),
    tui: bool = typer.Option(
        False, "--tui", help="Launch interactive terminal UI"
    ),
    multitarget: int | None = typer.Option(
        None,
        "--multitarget",
        help="Require candidate to be present in at least N input mRNAs "
        "(default: all input mRNAs)",
    ),
    mismatch: str = typer.Option(
        "1",
        "-m",
        "--mismatch",
        help="Stem mismatch position(s) on the passenger strand, comma/"
        "space-separated. Used both for scoring and, when oligos are "
        "designed, for the oligo's passenger arm. Pass 0 alone to disable "
        "mismatch introduction (perfect-complement passenger). Left unset, "
        "a named --backbone's own stored mismatch position(s) apply instead "
        "of this flag's default.",
    ),
    rna: bool = typer.Option(
        False, "--rna", help="Output sequences as RNA (U instead of T)"
    ),
    seed: int | None = typer.Option(
        None,
        "--seed",
        help="RNG seed for reproducibility (currently unused; reserved for "
        "future stochastic features). Echoed into metadata exports.",
    ),
    guide: str | None = typer.Option(
        None,
        "-g",
        "--guide",
        help="Design oligos for one guide (antisense) sequence directly, "
        "without scoring a FASTA.",
    ),
    backbone: str | None = typer.Option(
        None,
        "-b",
        "--backbone",
        help="Backbone name from the database. Given alongside INPUT, oligos "
        "are auto-designed for the resulting candidates. Required for "
        "--guide/--batch unless you define the insert yourself with --loop, "
        "--adapter-5 and --adapter-3.",
    ),
    batch: str | None = typer.Option(
        None,
        "--batch",
        help="CSV file of guide sequences to design oligos for directly, "
        "without scoring a FASTA (reads column 'antisense', 'guide', or "
        "'Guide Sequence').",
    ),
    loop: str | None = typer.Option(None, "--loop", help="Override loop sequence"),
    adapter_5: str | None = typer.Option(
        None, "--adapter-5", "--flank-5", help="Override 5' flank sequence"
    ),
    adapter_3: str | None = typer.Option(
        None, "--adapter-3", "--flank-3", help="Override 3' flank sequence"
    ),
    site_5: str | None = typer.Option(
        None,
        "--site-5",
        help=(
            "Override the backbone's left end: cut site 'G/AATTC' or its "
            "overhang 5'->3' — \"5'AATT\"/\"3'TGCA\"/'blunt'"
        ),
    ),
    site_3: str | None = typer.Option(
        None,
        "--site-3",
        help=(
            "Override the backbone's right end: cut site 'C/TCGAG' or its "
            "overhang 5'->3' — \"5'TCGA\"/\"3'GTAC\"/'blunt'"
        ),
    ),
    guide_on_5p: bool = typer.Option(
        False, "--guide-on-5p/--guide-on-3p", help="Place guide on 5' arm"
    ),
    keep_sites: bool | None = typer.Option(
        None,
        "--keep-sites/--no-keep-sites",
        help=(
            "Carry the rest of each recognition sequence on the insert, so "
            "ligation rebuilds the restriction site instead of scarring it "
            "(needs cut-site notation, not a bare overhang). Defaults to the "
            "backbone's own setting"
        ),
    ),
    list_backbones: bool = typer.Option(
        False, "--list-backbones", help="List available backbones and exit"
    ),
    restore_backbones: bool = typer.Option(
        False,
        "--restore-backbones",
        help="Bring back any bundled backbone you've deleted (e.g. from the "
        "TUI's 'd' key), leaving your own saved backbones untouched, and exit",
    ),
    info: bool = typer.Option(
        False,
        "--info",
        help="Show the version and a summary of the bundled reference data "
        "(coefficients, ΔG table, dG bins), and exit",
    ),
    version: bool = typer.Option(
        False,
        "--version",
        help="Show the version and exit",
    ),
) -> None:
    if version:
        console.print(f"KnockKnock {__version__}")
        return

    if info:
        _print_reference_info()
        return

    if list_backbones:
        _print_backbones_table()
        return

    if restore_backbones:
        restored = restore_bundled_backbones()
        if restored:
            console.print(f"[green]Restored:[/] {', '.join(restored)}")
        else:
            console.print("[dim]Nothing to restore — every bundled backbone is present.[/dim]")
        return

    oligo_only = guide is not None or batch is not None
    if guide is not None and batch is not None:
        console.print("[red]Error:[/] Use either --guide or --batch, not both.")
        raise typer.Exit(1)

    if oligo_only and fasta_inputs:
        console.print(
            "[red]Error:[/] INPUT (a FASTA to score) can't be combined with "
            "--guide/--batch, which design oligos directly from a guide you "
            "already have. Drop one or the other."
        )
        raise typer.Exit(1)

    if oligo_only:
        fmt = output_format or "fasta"
        if fmt not in ("fasta", "json"):
            console.print("[red]Error:[/] Unsupported format. Use fasta or json.")
            raise typer.Exit(1)

        config = _resolve_oligo_config(
            backbone,
            loop=loop,
            adapter_5=adapter_5,
            adapter_3=adapter_3,
            site_5=site_5,
            site_3=site_3,
            guide_on_5p=guide_on_5p,
            mismatch_positions=[1],
            keep_sites=keep_sites,
        )

        if batch:
            batch_path = Path(batch)
            if not batch_path.is_file():
                console.print(f"[red]Error:[/] File '{batch}' not found.")
                raise typer.Exit(1)

            try:
                batch_items = read_batch_csv(
                    batch_path.read_text(encoding="utf-8-sig"),
                    normalise_guide=lambda seq, label: _normalise_guide(
                        seq, label=label
                    ),
                )
            except (BatchInputError, ValueError) as exc:
                console.print(f"[red]Error:[/] {exc}")
                raise typer.Exit(1) from exc
            except OSError as exc:
                console.print(f"[red]Error:[/] Could not read '{batch}': {exc}")
                raise typer.Exit(1) from exc

            # A batch may contain user-supplied guides of different lengths.
            # Check against the shortest one, which guarantees every requested
            # mismatch position exists in every oligo we are about to assemble.
            try:
                config.mismatch_positions = _resolve_mismatch_positions(
                    ctx,
                    mismatch,
                    min(len(item.antisense) for item in batch_items),
                    backbone,
                )
            except ValueError as exc:
                console.print(f"[red]Error:[/] {exc}")
                raise typer.Exit(1) from exc

            pairs = assemble_batch(batch_items, config)

        else:
            try:
                clean_guide = _normalise_guide(guide, label="Guide")
                config.mismatch_positions = _resolve_mismatch_positions(
                    ctx, mismatch, len(clean_guide), backbone
                )
            except ValueError as exc:
                console.print(f"[red]Error:[/] {exc}")
                raise typer.Exit(1) from exc
            pairs = [assemble_oligos(config, clean_guide)]

        _warn_about_design(config, pairs)
        _emit_oligos(pairs, config, output_format=fmt, output=output, rna=rna)
        return

    # --- scoring mode ---
    inputs = fasta_inputs or ["-"]
    multi_input = len(inputs) > 1

    if inputs.count("-") > 1:
        console.print(
            "[red]Error:[/] '-' (stdin) can only be given once among multiple inputs."
        )
        raise typer.Exit(1)
    if multi_input and tui:
        console.print(
            "[red]Error:[/] --tui can only be used with a single input file."
        )
        raise typer.Exit(1)
    if multi_input and output not in (None, "-"):
        console.print(
            "[red]Error:[/] -o/--output cannot be used with multiple input "
            "files; each input gets its own auto-named output file. Omit "
            "-o, or pass a single input."
        )
        raise typer.Exit(1)

    if output_format is not None and output_format not in (
        "table",
        "fasta",
        "csv",
        "xlsx",
        "json",
    ):
        console.print(
            f"[red]Error:[/] Unsupported format '{output_format}'. "
            f"Use table, fasta, csv, xlsx, or json."
        )
        raise typer.Exit(1)

    # No explicit -f: pick based on whether stdout is a terminal. Interactive
    # use gets the rich table; piped/redirected use gets FASTA, since that's
    # what downstream sequence tools (BLAST, primer3, aligners) expect.
    # An explicit -f always wins, in either direction — e.g. -f table still
    # works when piped, -f csv still works at a terminal.
    resolved_format = output_format or ("table" if sys.stdout.isatty() else "fasta")

    # --backbone or any manual insert flag turns scoring into scoring+design:
    # the same candidates that get shown/exported also get cloning oligos.
    oligo_requested = bool(
        backbone
        or loop is not None
        or adapter_5 is not None
        or adapter_3 is not None
        or site_5 is not None
        or site_3 is not None
    )

    try:
        stem_mismatch = _resolve_mismatch_positions(ctx, mismatch, length, backbone)
    except ValueError as exc:
        console.print(f"[red]Error:[/] {exc}")
        raise typer.Exit(1) from exc

    # Resolved once, up front — the backbone/insert don't vary per input file
    # in a batch run, and failing fast here beats scoring every file first.
    oligo_config: OligoConfig | None = None
    if oligo_requested:
        oligo_config = _resolve_oligo_config(
            backbone,
            loop=loop,
            adapter_5=adapter_5,
            adapter_3=adapter_3,
            site_5=site_5,
            site_3=site_3,
            guide_on_5p=guide_on_5p,
            mismatch_positions=stem_mismatch,
            keep_sites=keep_sites,
        )

    # See --top's help: unlimited once oligo design is triggered, unless the
    # user typed -t explicitly (which always wins either way).
    effective_top_n = top_n if top_n is not None else (None if oligo_requested else 10)

    # Auto-named export paths are timestamped to the second (see
    # `export.resolve_export_path`), so two inputs sharing a stem (or the
    # same file listed twice) processed within the same second would
    # otherwise resolve to the same path and silently overwrite each other.
    # Tracked here, across the whole batch, so `_run_one` can disambiguate.
    seen_export_paths: set[Path] = set()

    for fasta_input in inputs:
        if multi_input:
            # stderr, so it never lands in a piped/redirected data stream —
            # just a progress marker for the human watching.
            err_console.print(f"[bold cyan]▶ {fasta_input}[/bold cyan]")
        _run_one(
            ctx,
            fasta_input,
            output=output,
            resolved_format=resolved_format,
            multi_input=multi_input,
            seen_export_paths=seen_export_paths,
            quality=quality,
            length=length,
            gc_min=gc_min,
            gc_max=gc_max,
            min_efficacy_score=min_efficacy_score,
            sort_dg=sort_dg,
            remove_motifs=remove_motifs,
            top_n=effective_top_n,
            tui=tui,
            multitarget=multitarget,
            mismatch=stem_mismatch,
            rna=rna,
            seed=seed,
            oligo_config=oligo_config,
        )


def _run_one(
    ctx: typer.Context,
    fasta_input: str,
    *,
    output: str | None,
    resolved_format: str,
    multi_input: bool,
    seen_export_paths: set[Path],
    quality: str,
    length: int,
    gc_min: float,
    gc_max: float,
    min_efficacy_score: float | None,
    sort_dg: bool,
    remove_motifs: bool,
    top_n: int | None,
    tui: bool,
    multitarget: int | None,
    mismatch: list[int],
    rna: bool,
    seed: int | None,
    oligo_config: OligoConfig | None,
) -> None:
    """Score, filter, and show/export/explore a single FASTA input.

    Split out of `run` so a multi-file invocation can call this once per
    input, unchanged, with each file's own `n_mrnas`-dependent defaults
    (`--multitarget`'s "all input mRNAs") resolved independently. When
    `oligo_config` is given, oligos are also designed for this file's own
    exported candidates.
    """
    mrna_tuples, input_warnings = _load_input(fasta_input)
    n_mrnas = len(mrna_tuples)

    # stderr, so a wildcard warning never lands in piped FASTA/CSV output.
    # --tui skips it: the TUI would wipe the line, and shows toasts instead.
    if not tui:
        for message in input_warnings:
            err_console.print(f"[yellow]Warning:[/] {message}")

    params = _resolve_params(
        ctx,
        n_mrnas=n_mrnas,
        length=length,
        quality=quality,
        gc_min=gc_min,
        gc_max=gc_max,
        min_efficacy_score=min_efficacy_score,
        sort_dg=sort_dg,
        remove_motifs=remove_motifs,
        top_n=top_n,
        multitarget=multitarget,
        mismatch=mismatch,
        rna=rna,
        seed=seed,
    )

    candidates_by_mrna = make_shrnas(mrna_tuples, params)
    results, stats = process_shrnas(candidates_by_mrna, mrna_tuples, params)

    if tui:
        tui_params = AnalysisParams(
            length=length,
            gc_min=0.0,
            gc_max=100.0,
            min_efficacy_score=None,
            remove_motifs=False,
            multitarget=1,
            rna_output=rna,
            stem_mismatch=mismatch,
        )
        tui_candidates, tui_stats = process_shrnas(
            candidates_by_mrna, mrna_tuples, tui_params
        )
        # The explorer starts from the unfiltered set, so an empty one here
        # means there is genuinely nothing to explore — open no UI at all
        # rather than an empty table with live filters that can't help.
        if not tui_candidates:
            for message in input_warnings:
                err_console.print(f"[yellow]Warning:[/] {message}")
            console.print(
                f"[yellow]{no_candidates_message(tui_stats, tui_params)}[/yellow]"
            )
            return
        tui_launch(
            tui_candidates,
            n_mrnas,
            fasta_input,
            rna_output=rna,
            warnings=input_warnings,
        )
        return

    # -t/--top caps every export format the same way it caps the table, so
    # copy-pasted (or piped/exported) candidates always match what -t asked
    # for and what the filters passed — not a separate "full dump" set. It
    # also caps how many candidates get oligos designed.
    exported = results[:top_n]

    pairs: list[OligoPair] = []
    if oligo_config is not None and exported:
        pairs = assemble_batch(exported, oligo_config)
        _warn_about_design(oligo_config, pairs)
    backbone_name = oligo_config.backbone_name if oligo_config is not None else None

    if resolved_format == "table":
        to_terminal(
            results,
            stats,
            params,
            top_n=top_n,
            total_mrnas=n_mrnas,
            input_path=fasta_input,
            oligo_pairs=pairs or None,
            backbone_name=backbone_name,
        )
        return

    # Non-table format: `-o -` (or no -o while piped) sends the data itself
    # to stdout and nothing else, so it's safe to feed to another program.
    # Batch mode (multiple inputs) ignores the tty check and always writes a
    # file per input instead — "one export per gene" is the whole point of
    # passing several files, and piping N files' worth of text into one
    # stdout stream isn't a usable substitute for that. `-o -` still forces
    # stdout explicitly if that's genuinely what's wanted.
    want_stdout = output == "-" or (
        output is None and not multi_input and not sys.stdout.isatty()
    )

    # An empty export is otherwise silent — an empty FASTA on a pipe, or a
    # CSV/XLSX with nothing but a header. Say why, on stderr so the data
    # stream stays clean. The `not want_stdout` branch below prints the table,
    # which carries the same message in its own panel.
    if not results and want_stdout:
        err_console.print(f"[yellow]{no_candidates_message(stats, params)}[/yellow]")

    if not want_stdout:
        # At a terminal with an explicit export requested: show the table
        # for a quick look, same as before, in addition to writing the file.
        to_terminal(
            results,
            stats,
            params,
            top_n=top_n,
            total_mrnas=n_mrnas,
            input_path=fasta_input,
            oligo_pairs=pairs or None,
            backbone_name=backbone_name,
        )

    ext = {"fasta": ".fasta", "csv": ".csv", "xlsx": ".xlsx", "json": ".json"}[
        resolved_format
    ]

    if resolved_format == "xlsx":
        if want_stdout:
            # A temp *directory*, not NamedTemporaryFile: xlsxwriter opens the
            # workbook by path, and Windows forbids a second open on a file
            # that still has an open handle. TemporaryDirectory has no such
            # restriction and cleans up the same way.
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir) / "candidates.xlsx"
                to_xlsx(
                    exported,
                    params,
                    tmp_path,
                    total_mrnas=n_mrnas,
                    input_path=fasta_input,
                    oligo_pairs=pairs or None,
                    backbone_name=backbone_name,
                )
                sys.stdout.buffer.write(tmp_path.read_bytes())
            sys.stdout.buffer.flush()
        else:
            out_path = _dedupe_export_path(
                _resolve_output_path(output, fasta_input, ext), seen_export_paths
            )
            to_xlsx(
                exported,
                params,
                out_path,
                total_mrnas=n_mrnas,
                input_path=fasta_input,
                oligo_pairs=pairs or None,
                backbone_name=backbone_name,
            )
            err_console.print(f"[green]Exported to {out_path}[/green]")
        return

    if resolved_format == "fasta":
        text = to_fasta(
            exported, params, oligo_pairs=pairs or None, backbone_name=backbone_name
        )
    elif resolved_format == "csv":
        text = to_csv_text(exported, params, input_path=fasta_input, total_mrnas=n_mrnas)
    elif resolved_format == "json":
        text = to_json_text(
            exported,
            params,
            input_path=fasta_input,
            oligo_pairs=pairs or None,
            backbone_name=backbone_name,
        )

    if want_stdout:
        sys.stdout.write(text)
        if pairs and resolved_format == "csv":
            err_console.print(
                "[yellow]Warning:[/] Oligo pairs were designed but CSV can't "
                "carry two tables over stdout; use -o to also write an "
                "accompanying *_oligos.csv file."
            )
    else:
        out_path = _dedupe_export_path(
            _resolve_output_path(output, fasta_input, ext), seen_export_paths
        )
        out_path.write_text(text, encoding="utf-8")
        err_console.print(f"[green]Exported to {out_path}[/green]")
        if pairs and resolved_format == "csv":
            oligo_path = _dedupe_export_path(
                _oligo_csv_path(output, fasta_input), seen_export_paths
            )
            oligo_path.write_text(
                export_oligos_csv_text(pairs, backbone_name or "Custom", rna=rna),
                encoding="utf-8",
            )
            err_console.print(
                f"[green]Exported {len(pairs)} oligo pair(s) to {oligo_path}[/green]"
            )


def _print_reference_info() -> None:
    """Show information about the bundled reference data."""
    coeffs = load_coeffs()
    deltaG = load_deltaG()
    dG_bins = load_dG_bins()

    console.print(f"[bold]KnockKnock[/bold] {__version__}\n")

    table = Table(title="Bundled Reference Data", expand=False)
    table.add_column("Data File", style="cyan")
    table.add_column("Details")

    coeffs_min = coeffs[["A", "G", "C", "T"]].min().min()
    coeffs_max = coeffs[["A", "G", "C", "T"]].max().max()
    table.add_row(
        "coeffs.csv",
        f"Shape: {coeffs.shape[0]} positions × {4} nucleotides\n"
        f"Value range: {coeffs_min:.3f} to {coeffs_max:.3f}\n"
        f"Source: Matveeva et al., 2012",
    )

    dg_vals = list(deltaG.values())
    table.add_row(
        "deltaG.csv",
        f"Entries: {len(deltaG)} dinucleotide pairs\n"
        f"Value range: {min(dg_vals):.2f} to {max(dg_vals):.2f}\n"
        f"Source: Xia et al., 1998",
    )

    dg_bins_min = dG_bins.min().min()
    dg_bins_max = dG_bins.max().max()
    table.add_row(
        "dG_bins.csv",
        f"Shape: {dG_bins.shape[0]} positions × {dG_bins.shape[1]} stability bins\n"
        f"Value range: {dg_bins_min:.2f} to {dg_bins_max:.2f}",
    )

    console.print(table)

    active_path = backbones_data_path()
    if active_path == bundled_backbones_path():
        console.print(
            f"\n[bold]Backbones:[/bold] {active_path} (bundled defaults — "
            "no custom edits yet)"
        )
    else:
        console.print(f"\n[bold]Backbones:[/bold] {active_path} (custom)")
