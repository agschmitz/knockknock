"""OligoScreen — backbone config + preview + batch export within the TUI."""

import shutil
import subprocess
from functools import lru_cache

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal
from textual.screen import Screen
from textual.widgets import Button, Checkbox, Header, Input, Label, Select, Static

from knockknock.backbone import (
    Backbone,
    backbone_label,
    delete_backbone,
    load_backbones,
    save_backbone,
)
from knockknock.cutsite import (
    BLUNT,
    LEFT,
    RIGHT,
    CutSiteError,
    Overhang,
    parse_overhang,
    validate_dna_sequence,
)
from knockknock.duplex_view import C_MM, render_diagram
from knockknock.engine import MAX_GUIDE_LENGTH, AnalysisParams, Candidate
from knockknock.export import resolve_export_path, write_candidates_export
from knockknock.export_modal import OLIGO_FORMATS, ExportFormatModal
from knockknock.help_screen import HelpScreen
from knockknock.modals import ConfirmModal, PromptModal
from knockknock.oligo import (
    OligoConfig,
    OligoPair,
    assemble_batch,
    assemble_oligos,
    compute_mismatched_passenger,
    export_oligos_csv_text,
    export_oligos_fasta_text,
    export_oligos_json_text,
    parse_mismatch_positions,
    validate_design,
)
from knockknock.two_row_footer import TwoRowFooter


@lru_cache(maxsize=1024)
def _cached_passenger(
    antisense: str, mismatch_positions: tuple[int, ...], guide_on_5p: bool
) -> str:
    """Memoize passenger generation so non-mismatch keystrokes don't recompute."""
    return compute_mismatched_passenger(antisense, list(mismatch_positions), guide_on_5p)


_CLIPBOARD_COMMANDS = (
    ("pbcopy",),
    ("wl-copy",),
    ("xclip", "-selection", "clipboard"),
    ("xsel", "--clipboard", "--input"),
    ("clip.exe",),
)


def _copy_to_system_clipboard(text: str) -> None:
    """Best-effort copy via whatever local clipboard tool is on PATH.

    Tried in addition to OSC 52 (see `OligoScreen.action_copy`) since that
    escape is a no-op on macOS Terminal.app. Every failure mode here —
    missing binary, no display, sandboxed/headless run — is expected and
    silently ignored; OSC 52 remains the fallback for those.
    """
    for command in _CLIPBOARD_COMMANDS:
        if shutil.which(command[0]) is None:
            continue
        try:
            subprocess.run(command, input=text.encode("utf-8"), check=True, timeout=2)
        except (OSError, subprocess.SubprocessError):
            continue
        return


_BLUNT_OVERHANG = Overhang(kind=BLUNT)


def _sites_kept(cfg: OligoConfig) -> str:
    """Summary fragment naming the sites `keep_sites` actually rebuilt.

    The diagram colours these bases but cannot label them — a one-base
    segment has no room for a name (see `duplex_view._render_label_line`) —
    so the summary line under the diagram is where they get spelled out.
    Empty when nothing was restored, so the line is unchanged by default.
    """
    restored = [
        f"{ovh.recognition} (+{cfg.site_bases(ovh)})"
        for ovh in (cfg.five_prime_overhang, cfg.three_prime_overhang)
        if cfg.site_bases(ovh)
    ]
    return f"  |  Sites kept: {', '.join(restored)}" if restored else ""


class SaveBackboneModal(PromptModal[str | None]):
    """Modal prompt for the name under which to save the current junctions/
    flanks/loop as a new row in the backbone database CSV."""

    DEFAULT_CSS = """
    SaveBackboneModal {
        align: center middle;
    }
    SaveBackboneModal > Container {
        width: auto;
        min-width: 44;
        height: auto;
        background: $surface;
        border: thick $primary;
        padding: 1 2;
    }
    SaveBackboneModal Input {
        width: 100%;
    }
    """

    def compose(self) -> ComposeResult:
        with Container():
            yield Static("Save the current backbone as:")
            yield Input(placeholder="Backbone name", id="backbone-name-input")
            yield Static("[dim]Enter to save, Esc to cancel[/dim]")

    def on_mount(self) -> None:
        self.query_one("#backbone-name-input", Input).focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.dismiss(event.value.strip() or None)


class OligoScreen(Screen):
    CSS = """
    #title {
        height: auto;
        padding: 1 2;
        text-style: bold;
        background: $surface;
        border-bottom: solid $primary;
    }

    #config-container {
        height: auto;
        padding: 1 2;
        border-bottom: solid $primary-background;
    }

    #config-container Horizontal {
        height: auto;
        margin: 1 0;
        align: center top;
    }

    #config-container Label {
        width: 14;
        padding: 0 1 0 0;
    }

    #config-container Input {
        width: 24;
    }

    #config-container Select {
        width: 34;
    }

    #keep-sites {
        margin: 1 0 0 2;
    }

    #save-backbone-btn {
        /* Select wraps its visible control in a 1-row margin on both sides
        (for the dropdown overlay), so nudge the button down to match its
        actual box rather than the Select widget's outer (taller) region. */
        margin: 1 0 0 2;
        min-width: 18;
    }

    #save-backbone-btn.-hidden {
        display: none;
    }

    #preview-container {
        height: 1fr;
        padding: 1 2;
        border-bottom: solid $primary-background;
        overflow-y: auto;
    }

    #preview-text {
        height: auto;
        content-align-horizontal: center;
        width: 100%;
    }
    """

    BINDINGS = [
        # Design actions (row 1 of the footer — see `row_break` in compose()).
        Binding("e", "export", "Export All"),
        Binding("s", "save_backbone", "Save Backbone"),
        Binding("d", "delete_backbone", "Delete Backbone"),
        Binding("c", "copy", "Copy Example"),
        # Navigation / app (row 2).
        Binding("left", "prev_cand", "Prev"),
        Binding("right", "next_cand", "Next"),
        Binding("b", "go_back", "Back"),
        Binding("question_mark", "help", "Help"),
        Binding("q", "quit", "Quit"),
        Binding("escape", "unfocus", "Unfocus", show=False),
    ]

    def __init__(
        self,
        candidates: list[Candidate],
        backbone_name: str | None = None,
        input_path: str = "",
        params: AnalysisParams | None = None,
        total_mrnas: int = 1,
    ) -> None:
        super().__init__()
        self._candidates = candidates
        self._candidate_idx = 0
        self._input_path = input_path
        self._params = params or AnalysisParams()
        self._total_mrnas = total_mrnas
        self._backbones = load_backbones()
        self._config_errors: list[str] = []
        self._selected_backbone_name = backbone_name or (
            self._backbones[0].name if self._backbones else "Custom"
        )

        bb = next(
            (b for b in self._backbones if b.name == self._selected_backbone_name),
            None,
        )
        if bb is not None:
            self._config = OligoConfig.from_backbone(bb)
        elif self._backbones:
            bb = self._backbones[0]
            self._config = OligoConfig.from_backbone(bb)
            self._selected_backbone_name = bb.name
        else:
            self._config = OligoConfig(
                five_prime_overhang=_BLUNT_OVERHANG,
                five_prime_flank="",
                loop="",
                three_prime_flank="",
                three_prime_overhang=_BLUNT_OVERHANG,
            )
            self._selected_backbone_name = "Custom"

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)

        with Container(id="title"):
            yield Static(
                f"Oligo Designer  —  {len(self._candidates)} candidates selected"
            )

        with Container(id="config-container"):
            # Backbone selector
            with Horizontal():
                yield Label("Backbone:")
                opts = [(backbone_label(b), b.name) for b in self._backbones] + [
                    ("Custom", "Custom")
                ]
                yield Select(options=opts, id="backbone-select")
                yield Button("Save Backbone", id="save-backbone-btn", variant="primary")

            # 5' row
            with Horizontal():
                yield Label("5' Site:")
                yield Input(placeholder="e.g. G/AATTC", id="site-5")
                yield Label("5' Flank:")
                yield Input(id="flank-5")

            # Loop + mismatch row
            with Horizontal():
                yield Label("Loop:")
                yield Input(id="loop")
                yield Label("Mismatch pos:")
                yield Input(id="mismatch")

            # 3' row
            with Horizontal():
                yield Label("3' Site:")
                yield Input(placeholder="e.g. C/TCGAG", id="site-3")
                yield Label("3' Flank:")
                yield Input(id="flank-3")

            # Orientation
            with Horizontal():
                yield Label("Guide on:")
                yield Select(
                    options=[("3' arm", "3p"), ("5' arm", "5p")],
                    id="guide-arm",
                )
                yield Checkbox("Keep restriction sites", id="keep-sites")

        with Container(id="preview-container"):
            yield Static("", id="preview-text")

        yield TwoRowFooter(row_break="prev_cand")

    def on_mount(self) -> None:
        # Programmatic sync, not a user edit — must not trip _mark_custom()
        # via the Input.Changed / Select.Changed handlers below.
        with self.prevent(Input.Changed, Select.Changed, Checkbox.Changed):
            self.query_one(
                "#backbone-select", Select
            ).value = self._selected_backbone_name
            self.query_one("#site-5", Input).value = self._config.five_prime_overhang.spec
            self.query_one("#flank-5", Input).value = self._config.five_prime_flank
            self.query_one("#loop", Input).value = self._config.loop
            self.query_one("#mismatch", Input).value = ",".join(
                str(m) for m in self._config.mismatch_positions
            )
            self.query_one(
                "#site-3", Input
            ).value = self._config.three_prime_overhang.spec
            self.query_one("#flank-3", Input).value = self._config.three_prime_flank
            guide_arm_val = "5p" if self._config.guide_on_5p else "3p"
            self.query_one("#guide-arm", Select).value = guide_arm_val
            self.query_one("#keep-sites", Checkbox).value = self._config.keep_sites
        self._update_backbone_button()
        self._refresh_preview()

    def _update_backbone_button(self) -> None:
        """One button, two jobs, decided by what is selected.

        "Custom" means the fields have diverged from any saved backbone, so the
        useful action is to save them. A named backbone is selected means the
        useful action is to remove it. The button is never hidden, so there is
        always a visible, mouse-reachable way to manage the database.
        """
        button = self.query_one("#save-backbone-btn", Button)
        if self._selected_backbone_name == "Custom":
            button.label = "Save Backbone"
            button.variant = "primary"
        else:
            button.label = "Delete Backbone"
            button.variant = "error"

    # ------------------------------------------------------------------
    #  Config helpers
    # ------------------------------------------------------------------

    def _read_site(self, widget_id: str, label: str) -> Overhang:
        """Parse a junction input, recording (not raising) any parse error.

        ``junction`` is not cosmetic: a directly-stated overhang describes the
        backbone's own end, and the backbone protrudes with opposite strands
        at the two ends — so the same ``5'TTAT`` means different things in the
        two fields (see `cutsite.parse_overhang`).
        """
        raw = self.query_one(f"#{widget_id}", Input).value
        junction = LEFT if widget_id == "site-5" else RIGHT
        try:
            return parse_overhang(raw, junction=junction)
        except CutSiteError as exc:
            self._config_errors.append(f"{label}: {exc}")
            return Overhang(kind=BLUNT, junction=junction)

    def _read_sequence(self, widget_id: str, label: str) -> str:
        """Read a free-text DNA field, recording (not raising) any bad base."""
        raw = self.query_one(f"#{widget_id}", Input).value.strip().upper()
        try:
            validate_dna_sequence(raw, label)
        except CutSiteError as exc:
            self._config_errors.append(str(exc))
            return ""
        return raw

    def _selected_backbone_label(self) -> str:
        """Display name for the currently-selected backbone, enzymes appended.

        Falls back to the raw selected name ("Custom", or a name that no
        longer matches a loaded backbone) when there is nothing to look up.
        """
        bb = next(
            (b for b in self._backbones if b.name == self._selected_backbone_name),
            None,
        )
        return backbone_label(bb) if bb is not None else self._selected_backbone_name

    def _read_config(self) -> OligoConfig:
        self._config_errors = []
        left = self._read_site("site-5", "5' site")
        right = self._read_site("site-3", "3' site")
        flank5 = self._read_sequence("flank-5", "5' flank")
        loop = self._read_sequence("loop", "Loop")
        flank3 = self._read_sequence("flank-3", "3' flank")

        # Same strict parser the CLI uses. This screen used to hand-parse and
        # silently drop anything unparseable, so "1,99,abc" was a hard error on
        # the command line and a silent [1] here. Parse errors now land in
        # _config_errors and show under the diagram, like every other bad field.
        raw = self.query_one("#mismatch", Input).value.strip()
        guide = self._current_candidate()
        guide_length = len(guide.antisense) if guide is not None else MAX_GUIDE_LENGTH
        try:
            mismatches = parse_mismatch_positions(raw, guide_length)
        except ValueError as exc:
            self._config_errors.append(f"Mismatch positions {exc}")
            mismatches = [1]

        guide_on_5p = self.query_one("#guide-arm", Select).value == "5p"
        keep_sites = self.query_one("#keep-sites", Checkbox).value

        return OligoConfig(
            five_prime_overhang=left,
            five_prime_flank=flank5,
            loop=loop,
            three_prime_flank=flank3,
            three_prime_overhang=right,
            guide_on_5p=guide_on_5p,
            mismatch_positions=mismatches,
            backbone_name=self._selected_backbone_label(),
            keep_sites=keep_sites,
        )

    def _fill_from_backbone(self, bb) -> None:
        # Programmatic sync, not a user edit — must not trip _mark_custom()
        # via the Input.Changed / Select.Changed handlers below (which would
        # otherwise immediately flip the freshly-picked backbone to "Custom").
        with self.prevent(Input.Changed, Select.Changed, Checkbox.Changed):
            self.query_one("#site-5", Input).value = bb.five_prime_overhang.spec
            self.query_one("#flank-5", Input).value = bb.five_prime_flank
            self.query_one("#loop", Input).value = bb.loop
            self.query_one("#mismatch", Input).value = ",".join(
                str(m) for m in bb.mismatch_positions
            )
            self.query_one("#site-3", Input).value = bb.three_prime_overhang.spec
            self.query_one("#flank-3", Input).value = bb.three_prime_flank
            self.query_one("#guide-arm", Select).value = "5p" if bb.guide_on_5p else "3p"
            self.query_one("#keep-sites", Checkbox).value = bb.keep_sites

    # ------------------------------------------------------------------
    #  Event handlers
    # ------------------------------------------------------------------

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "backbone-select":
            name = event.value
            if name == "Custom":
                self._selected_backbone_name = "Custom"
            else:
                bb = next((b for b in self._backbones if b.name == name), None)
                if bb is not None:
                    self._selected_backbone_name = bb.name
                    self._fill_from_backbone(bb)
            self._update_backbone_button()
            self._config = self._read_config()
            self._refresh_preview()
        elif event.select.id == "guide-arm":
            self._mark_custom()
            self._config = self._read_config()
            self._refresh_preview()

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        if event.checkbox.id == "keep-sites":
            # `keep_sites` is a stored backbone column, so changing it really
            # does mean "this is no longer that backbone" — same as editing a
            # site or flank field. Programmatic syncs are shielded by the
            # `prevent(..., Checkbox.Changed)` blocks, so this only ever fires
            # for a real click or keypress.
            self._mark_custom()
            self._config = self._read_config()
            self._refresh_preview()

    def on_input_changed(self, event: Input.Changed) -> None:
        if event.input.id in (
            "site-5",
            "flank-5",
            "loop",
            "site-3",
            "flank-3",
            "mismatch",
        ):
            self._mark_custom()
            self._config = self._read_config()
            self._refresh_preview()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id != "save-backbone-btn":
            return
        if self._selected_backbone_name == "Custom":
            self.action_save_backbone()
        else:
            self.action_delete_backbone()

    def _mark_custom(self) -> None:
        if self._selected_backbone_name != "Custom":
            self._selected_backbone_name = "Custom"
            self.query_one("#backbone-select", Select).value = "Custom"
            self._update_backbone_button()

    # ------------------------------------------------------------------
    #  Preview rendering
    # ------------------------------------------------------------------

    # The diagram itself lives in `knockknock.duplex_view` — it is a pure
    # function of the config and the two arms, so it stays testable (and
    # reusable by a non-terminal front end) without this screen. What is left
    # here is only what needs widget state: reading the fields, picking the
    # current candidate, and pushing the rendered rows into the Static.

    def _current_candidate(self) -> Candidate | None:
        if not self._candidates:
            return None
        return self._candidates[self._candidate_idx]

    def _current_pair(self) -> tuple[Candidate, OligoConfig, OligoPair] | None:
        """Re-read the config fields and assemble the pair for the candidate
        on screen. Shared by the preview and the clipboard-copy action so
        they can't drift on how the current pair is derived."""
        cand = self._current_candidate()
        if cand is None:
            return None
        self._config = self._read_config()
        passenger = _cached_passenger(
            cand.antisense,
            tuple(self._config.mismatch_positions),
            self._config.guide_on_5p,
        )
        pair = assemble_oligos(self._config, cand.antisense, passenger, name=cand.name)
        return cand, self._config, pair

    def _refresh_preview(self) -> None:
        preview = self.query_one("#preview-text", Static)

        current = self._current_pair()
        if current is None:
            preview.update("No candidates selected.")
            return
        cand, cfg, pair = current
        passenger = pair.passenger

        arm1, arm2 = cfg.arms(cand.antisense, passenger)

        row1, row2, row3, row4 = render_diagram(cfg, arm1, arm2)

        lines: list[str] = [
            f"[bold]Candidate {self._candidate_idx + 1}/{len(self._candidates)}: "
            f"{cand.name}[/bold]\n",
            f"Guide:     {cand.antisense}",
            f"Passenger: {passenger}\n",
            row1,
            row2,
            row3,
            f"[dim]{row4}[/dim]\n",
            f"Backbone ends: {cfg.five_prime_overhang.describe()} / "
            f"{cfg.three_prime_overhang.describe()}{_sites_kept(cfg)}  |  "
            f"Annealing Tm: {pair.tm:.1f}°C  |  Overlap: {pair.overlap} bp  |  "
            f"Fwd: {len(pair.forward)} nt  Rev: {len(pair.reverse)} nt",
        ]

        for message in self._config_errors:
            lines.append(f"[{C_MM}]✗ {message}[/{C_MM}]")
        for message in validate_design(cfg, arm1, arm2):
            lines.append(f"[yellow]⚠ {message}[/yellow]")

        preview.update("\n".join(lines))

    # ------------------------------------------------------------------
    #  Actions
    # ------------------------------------------------------------------

    def action_export(self) -> None:
        if not self._candidates:
            self.notify("No candidates to export.", severity="warning")
            return
        self.app.push_screen(
            ExportFormatModal(OLIGO_FORMATS, title="Export candidates + oligo pairs as:"),
            callback=self._on_export_format_chosen,
        )

    def _on_export_format_chosen(self, fmt: str | None) -> None:
        """Write the candidates that went into this batch to their own file,
        alongside a second file with the assembled oligo pairs.

        The candidates file reuses `export.write_candidates_export` — the
        same dispatch `tui.py`'s explorer screen calls for its own `e`
        export — so the two screens can't drift on what a "candidates
        export" contains for a given format. The oligo pairs file reuses
        this module's own `export_oligos_*_text` formatters, unchanged from
        before. Two files, every format, rather than one combined file —
        simpler than threading an `oligo_pairs` section through the
        candidate exporters, and it means a CSV, FASTA, or JSON candidates
        file exported here is byte-identical in shape to one exported from
        the explorer screen.
        """
        if fmt is None:
            return

        self._config = self._read_config()
        pairs = assemble_batch(self._candidates, self._config)
        name = self._config.backbone_name

        ext = {f.value: f.ext for f in OLIGO_FORMATS}[fmt]
        candidates_path = resolve_export_path(
            ext, reference_path=self._input_path, stem="oligos", label="candidates"
        )
        oligos_path = resolve_export_path(
            ext, reference_path=self._input_path, stem="oligos", label="oligos"
        )

        if fmt == "csv":
            oligos_text = export_oligos_csv_text(pairs, name)
        elif fmt == "json":
            oligos_text = export_oligos_json_text(pairs, name)
        else:  # fasta
            oligos_text = export_oligos_fasta_text(pairs, name)

        def _write() -> None:
            write_candidates_export(
                fmt,
                self._candidates,
                self._params,
                candidates_path,
                input_path=self._input_path,
                total_mrnas=self._total_mrnas,
            )
            oligos_path.write_text(oligos_text, encoding="utf-8")

        try:
            _write()
        except OSError as exc:
            self.notify(f"Could not write export: {exc}", severity="error")
            return
        self.notify(
            f"Exported {len(self._candidates)} candidates to "
            f"[bold]{candidates_path}[/bold] and {len(pairs)} oligo pairs to "
            f"[bold]{oligos_path}[/bold]",
            timeout=5,
        )

    def action_save_backbone(self) -> None:
        self._config = self._read_config()
        if self._config_errors:
            self.notify(
                "Fix the highlighted errors before saving this as a backbone.",
                severity="error",
            )
            return
        self.app.push_screen(SaveBackboneModal(), callback=self._on_backbone_named)

    def action_delete_backbone(self) -> None:
        """Remove the selected backbone, after a confirmation prompt."""
        name = self._selected_backbone_name
        if name == "Custom":
            self.notify(
                "Select a saved backbone to delete. 'Custom' is not stored.",
                severity="warning",
            )
            return
        self.app.push_screen(
            ConfirmModal(
                f"Delete the backbone '{name}'?",
                detail=(
                    "Bundled backbones can be brought back with "
                    "'knockknock --restore-backbones'."
                ),
            ),
            callback=lambda ok: self._on_delete_confirmed(name, ok),
        )

    def _on_delete_confirmed(self, name: str, confirmed: bool | None) -> None:
        if not confirmed:
            return
        try:
            delete_backbone(name)
        except ValueError as exc:
            self.notify(str(exc), severity="error")
            return

        self._backbones = load_backbones()
        self._selected_backbone_name = "Custom"
        opts = [(backbone_label(b), b.name) for b in self._backbones] + [("Custom", "Custom")]
        with self.prevent(Input.Changed, Select.Changed, Checkbox.Changed):
            select = self.query_one("#backbone-select", Select)
            select.set_options(opts)
            select.value = "Custom"
        self._update_backbone_button()
        self._refresh_preview()
        self.notify(f"Deleted backbone [bold]{name}[/bold].", timeout=5)

    def _on_backbone_named(self, name: str | None) -> None:
        if not name:
            return
        cfg = self._config
        backbone = Backbone(
            name=name,
            five_prime_enzyme="",
            five_prime_overhang=cfg.five_prime_overhang,
            five_prime_flank=cfg.five_prime_flank,
            loop=cfg.loop,
            three_prime_flank=cfg.three_prime_flank,
            three_prime_enzyme="",
            three_prime_overhang=cfg.three_prime_overhang,
            guide_on_5p=cfg.guide_on_5p,
            keep_sites=cfg.keep_sites,
            mismatch_positions=cfg.mismatch_positions,
        )
        try:
            save_backbone(backbone)
        except ValueError as exc:
            self.notify(str(exc), severity="error")
            return

        self._backbones = load_backbones()
        self._selected_backbone_name = name
        opts = [(backbone_label(b), b.name) for b in self._backbones] + [("Custom", "Custom")]
        with self.prevent(Input.Changed, Select.Changed, Checkbox.Changed):
            select = self.query_one("#backbone-select", Select)
            select.set_options(opts)
            select.value = name
        self._update_backbone_button()
        self.notify(
            f"Saved backbone [bold]{name}[/bold] to the backbone database.",
            timeout=5,
        )

    def action_copy(self) -> None:
        """Copy the current pair to the system clipboard.

        Text comes from `export_oligos_fasta_text` — the same formatter the
        FASTA export uses — so naming and layout can't drift between copy and
        export; this only supplies the single pair and the current backbone.

        Sends Textual's own OSC 52 escape (``app.copy_to_clipboard``) so a
        remote/SSH terminal that honours it picks the text up, *and* shells
        out to whatever local clipboard tool is on PATH. OSC 52 alone silently
        does nothing on macOS Terminal.app — Textual's own docs call this out
        — which was the everyday case this left broken; the shell-out is the
        fallback that actually reaches the clipboard there (and on plain X11).
        Both are attempted; either one succeeding is enough.
        """
        current = self._current_pair()
        if current is None:
            return
        _cand, cfg, pair = current
        text = export_oligos_fasta_text([pair], cfg.backbone_name)

        self.app.copy_to_clipboard(text)
        _copy_to_system_clipboard(text)
        self.notify("Copied example oligo pair to clipboard", timeout=2)

    def action_next_cand(self) -> None:
        if self._candidates:
            self._candidate_idx = (self._candidate_idx + 1) % len(self._candidates)
            self._refresh_preview()

    def action_prev_cand(self) -> None:
        if self._candidates:
            self._candidate_idx = (self._candidate_idx - 1) % len(self._candidates)
            self._refresh_preview()

    def action_go_back(self) -> None:
        self.app.pop_screen()

    def action_unfocus(self) -> None:
        """Escape out of a config field without leaving the screen.

        Lets single-letter shortcuts (`e`/`s`/`d`/`c`/`b`) work again right
        after editing a field, instead of popping back to the candidate
        table like `b` does.
        """
        if isinstance(self.focused, (Input, Select)):
            self.set_focus(None)

    def action_help(self) -> None:
        self.app.push_screen(HelpScreen())
