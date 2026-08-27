"""Textual interactive terminal UI (--tui)."""

import json
import shlex
import webbrowser
from collections.abc import Callable
from pathlib import Path
from typing import Literal
from urllib.parse import quote

from rich.text import Text
from textual import events, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal, Vertical
from textual.css.query import NoMatches
from textual.screen import Screen
from textual.timer import Timer
from textual.widgets import (
    Button,
    Checkbox,
    DataTable,
    DirectoryTree,
    Footer,
    Header,
    Input,
    LoadingIndicator,
    Static,
    TextArea,
    Tree,
)
from textual.widgets._data_table import CellDoesNotExist
from textual.widgets._directory_tree import DirEntry
from textual.widgets._tree import TreeNode

from knockknock.engine import (
    AnalysisParameterError,
    AnalysisParams,
    Candidate,
    candidate_sort_key,
    make_shrnas,
    process_shrnas,
    rejection_reason,
    validate_analysis_params,
)
from knockknock.export import (
    format_candidate_cells,
    no_candidates_message,
    resolve_export_path,
    write_candidates_export,
)
from knockknock.export_modal import (
    CANDIDATE_FORMATS,
    ExportFormatModal,
    export_and_notify,
)
from knockknock.help_screen import HelpScreen
from knockknock.modals import PromptModal
from knockknock.oligo_tui import OligoScreen
from knockknock.seqinput import (
    SequenceInputError,
    load_sequences,
    read_sequence_file,
)
from knockknock.two_row_footer import TwoRowFooter

_STATE_FILE = Path.home() / ".config" / "knockknock" / "state.json"

#: Below this many candidates passing the default filters (efficacy score
#: ≥0, GC 30–60%, all targets), `ExploreApp` relaxes them rather than opening
#: on a near-empty table.
_MIN_DEFAULT_RESULTS = 8

# NCBI's interactive BLAST form (as opposed to its CMD=Put/Get automated API)
# happens to read its own GET parameters back to prefill the query box and
# algorithm — the same trick NCBI's own "BLAST this sequence" links use. No
# DATABASE param on purpose: NCBI remembers the last database a given browser
# searched against (a cookie) and that saved choice silently overrides
# whatever DATABASE we'd put in the URL after the first-ever search anyway,
# so specifying one here would just be misleading. BLAST_PROGRAMS=blastn picks
# the "somewhat similar sequences" algorithm rather than the default
# Megablast, which is tuned for near-identical long sequences and a poor fit
# for ~19-25 nt guides. Sequence count is capped well short of any URL-length
# limit.
_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_BLAST_MAX_SEQS = 25


def _load_state() -> dict:
    """Read the UI state file, treating any problem as "no saved state".

    OSError covers a missing/unreadable file, ValueError covers corrupt JSON
    (json.JSONDecodeError subclasses it). Narrow on purpose: a TypeError or
    AttributeError here would be a bug in our own code and should surface.
    """
    try:
        return json.loads(_STATE_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _save_state(state: dict) -> None:
    try:
        _STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
        _STATE_FILE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    except OSError:
        pass  # best-effort: a read-only home must not crash the TUI


def _ctrl_key_display(binding: Binding) -> str | None:
    """`ctrl+o` -> `Ctrl+O`, instead of Textual's default `^o` caret.

    Returns `None` for anything else (an explicit `key_display` override,
    a non-ctrl key, a ctrl combo with more than one modifier/key) so the
    caller falls back to `App.get_key_display`'s own formatting.
    """
    if binding.key_display:
        return None
    modifiers, key = binding.parse_key()
    if modifiers == ["ctrl"] and len(key) == 1:
        return f"Ctrl+{key.upper()}"
    return None


class TopNModal(PromptModal[int | None]):
    """Modal input dialog for selecting the top N guides."""

    DEFAULT_CSS = """
    TopNModal {
        align: center middle;
    }
    TopNModal > Container {
        width: auto;
        min-width: 40;
        height: auto;
        background: $surface;
        border: thick $primary;
        padding: 1 2;
    }
    TopNModal Input {
        width: 100%;
    }
    """

    def compose(self) -> ComposeResult:
        with Container():
            yield Static("Select top N guides from the current view:")
            yield Input(placeholder="e.g. 4, 8, 10", id="top-n-input")
            yield Static("[dim]Enter to confirm, Esc to cancel[/dim]")

    def on_mount(self) -> None:
        self.query_one("#top-n-input", Input).focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self._submit_value(event.value.strip())

    def _submit_value(self, raw: str) -> None:
        try:
            n = int(raw)
            if n > 0:
                self.dismiss(n)
                return
        except ValueError:
            pass
        self.dismiss(None)


class ExploreApp(App):
    """Interactive shRNA candidate explorer."""

    TITLE = "KnockKnock"
    SUB_TITLE = "shRNA candidate explorer"

    CSS = """
    Screen {
        layout: vertical;
    }

    /* Textual's default border style for Input/Select/Checkbox ("tall")
    renders the side edges as partial-height block glyphs that many
    terminal fonts draw with visible gaps between rows. This App's CSS
    cascades to every screen pushed onto it (TopNModal, OligoScreen,
    SaveBackboneModal), so overriding it once here to a plain single-line
    border keeps every entry field, dropdown, and checkbox in the TUI a
    clean rectangle. */
    Input, SelectCurrent, SelectOverlay, ToggleButton {
        border: solid $primary-background;
    }

    Input:focus, Select:focus > SelectCurrent, ToggleButton:focus {
        border: solid $primary;
    }

    #onboarding-banner {
        height: auto;
        background: $accent 20%;
        color: $text;
        padding: 0 2;
    }

    #summary-bar {
        height: auto;
        background: $surface;
        padding: 1 2;
        border-bottom: solid $primary;
    }

    #summary-bar Static {
        width: 1fr;
    }

    #summary-counts {
        text-style: bold;
    }

    DataTable {
        height: 1fr;
        border: solid $primary-background;
    }

    #filter-bar {
        height: auto;
        background: $boost;
        padding: 1 2;
        border-top: solid $primary;
        align: center top;
    }

    .filter-group {
        width: auto;
        height: auto;
        margin: 0 2;
        align: center top;
    }

    .score-group {
        border: solid $accent;
        padding: 0 1;
    }

    .score-group Static {
        color: $accent;
        text-style: bold;
    }

    #filter-bar Input {
        width: 14;
    }

    #filter-bar Checkbox {
        width: auto;
    }

    #design-oligos-btn {
        margin-top: 1;
    }
    """

    BINDINGS = [
        Binding("t", "toggle_target_all", "Target all"),
        Binding("T", "focus_score", "Efficacy min"),
        Binding("g", "focus_gc_min", "GC min"),
        Binding("G", "focus_gc_max", "GC max"),
        Binding("m", "toggle_motifs", "Low complexity motifs"),
        Binding("s", "toggle_sort", "Sort mode"),
        Binding("space", "toggle_selection", "Select row"),
        Binding("a", "select_all_filtered", "Select all"),
        Binding("A", "clear_selection", "Clear"),
        Binding("n", "select_top_n", "Top N"),
        Binding("o", "design_oligos", "Design Oligos"),
        Binding("b", "blast_selected", "BLAST"),
        Binding("e", "export_view", "Export"),
        Binding("ctrl+o", "open_file", "Open file"),
        Binding("question_mark", "help", "Help"),
        Binding("x", "dismiss_banner", "Dismiss banner", show=False),
        Binding("q", "quit", "Quit"),
        Binding("escape", "unfocus", "Unfocus", show=False),
    ]

    def __init__(
        self,
        all_candidates: list[Candidate],
        total_mrnas: int,
        input_path: str,
        rna_output: bool = False,
        warnings: list[str] | None = None,
    ):
        super().__init__()
        self._all_candidates = all_candidates
        self._total_mrnas = total_mrnas
        self._input_path = input_path
        self._rna_output = rna_output
        #: Input notes (wildcard bases) raised while loading. Shown as toasts
        #: on mount rather than printed by the caller: the app takes over the
        #: screen immediately, so anything written before `run()` is wiped.
        self._warnings = warnings or []
        self._sort_by_dg = False
        self._motif_filter = True
        self._target_all_only = True
        self._gc_min = 30.0
        self._gc_max = 60.0
        self._score_min: float | None = 0.0
        self._filter_error: str | None = None
        self._selected: set[str] = set()  # candidate names checked for oligo design
        self._refresh_timer: Timer | None = (
            None  # debounce handle for filter-input changes
        )
        #: Set by ``ctrl+o`` before exiting, and read by `launch` — an
        #: `App` can't push another `App`, so "open a different file" has to
        #: be a request to the loop that owns both apps rather than something
        #: this screen does itself.
        self.open_requested = False

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        state = _load_state()
        if not state.get("seen_banner"):
            yield Static(
                "Welcome — filter at the bottom · [b]space[/b] select · "
                "[b]o[/b] design oligos · [b]e[/b] export · [b]?[/b] help · "
                "[b]q[/b] quit  [dim]([b]x[/b] dismiss)[/dim]",
                id="onboarding-banner",
            )
        with Container(id="summary-bar"):
            yield Static("", id="summary-input")
            yield Static("", id="summary-filters")
            yield Static("", id="summary-counts")
        yield DataTable(id="results-table", cursor_type="row")
        with Horizontal(id="filter-bar"):
            with Vertical(classes="filter-group score-group"):
                yield Static("Efficacy \u2265")
                yield Input(
                    value="0",
                    placeholder="any",
                    id="score-min",
                )
            with Vertical(classes="filter-group"):
                yield Static("GC Min")
                yield Input(
                    value="30",
                    placeholder="%",
                    id="gc-min",
                )
            with Vertical(classes="filter-group"):
                yield Static("GC Max")
                yield Input(
                    value="60",
                    placeholder="%",
                    id="gc-max",
                )
            with Vertical(classes="filter-group"):
                yield Static(" ")
                yield Checkbox(
                    "Filter low complexity motifs", value=True, id="motif-check"
                )
                yield Checkbox("Sort by ΔG Score", value=False, id="sort-check")
                yield Checkbox("Target all inputs", value=True, id="target-all-check")
            with Vertical(classes="filter-group"):
                yield Static(" ")
                yield Button(
                    "Design Oligos", id="design-oligos-btn", variant="primary"
                )
        # Row break lands between the filter/selection cluster and the
        # actions/app cluster in BINDINGS below — keep them in sync.
        yield TwoRowFooter(row_break="design_oligos")

    # Action names that should not fire while a text Input has focus,
    # so single-letter shortcuts don't conflict with typing in the filter bar.
    _INPUT_GATED_ACTIONS: frozenset[str] = frozenset(
        {
            "toggle_selection",
            "select_all_filtered",
            "clear_selection",
            "select_top_n",
            "toggle_motifs",
            "toggle_target_all",
            "toggle_sort",
            "design_oligos",
            "blast_selected",
            "export_view",
            "dismiss_banner",
        }
    )

    # These are bound on the App (not a Screen) so the single-letter shortcuts
    # keep working regardless of which Input has focus — but that also means
    # they stay "active" (and show up in the Footer) once a Screen like the
    # oligo designer is pushed on top, even though they only apply to the
    # candidate table. Gated back down to just the table view here.
    _TABLE_VIEW_ACTIONS: frozenset[str] = frozenset(
        {
            "focus_score",
            "focus_gc_min",
            "focus_gc_max",
            "toggle_motifs",
            "toggle_target_all",
            "toggle_sort",
            "toggle_selection",
            "select_all_filtered",
            "clear_selection",
            "select_top_n",
            "design_oligos",
            "blast_selected",
            "export_view",
            "dismiss_banner",
        }
    )

    def check_action(self, action: str, parameters: tuple) -> bool | None:
        if action in self._TABLE_VIEW_ACTIONS and isinstance(self.screen, OligoScreen):
            return False  # table-only shortcut; not shown/usable on this screen
        if action in self._INPUT_GATED_ACTIONS and isinstance(self.focused, Input):
            return False  # hidden from footer + binding suppressed
        return True

    def get_key_display(self, binding: Binding) -> str:
        return _ctrl_key_display(binding) or super().get_key_display(binding)

    def on_mount(self) -> None:
        table = self.query_one("#results-table", DataTable)
        table.add_columns(
            "Sel",
            "Rank",
            "Guide Sequence",
            "Efficacy Score",
            "ΔG Score",
            "GC%",
            "ΔG",
            "Targets",
            "Pos",
            "mRNA",
        )
        if len(self._filtered()) < _MIN_DEFAULT_RESULTS:
            self._relax_default_filters()
        self._refresh()
        for message in self._warnings:
            self.notify(message, title="Input warning", severity="warning", timeout=10)

    def _relax_default_filters(self) -> None:
        """Fall back to wide-open filters when the default ones (efficacy
        score ≥0, GC 30–60%, all targets) are too strict for this input to
        clear the bar of `_MIN_DEFAULT_RESULTS` candidates."""
        self._gc_min = 0.0
        self._gc_max = 100.0
        self._score_min = None
        self._target_all_only = False
        self.query_one("#gc-min", Input).value = "0"
        self.query_one("#gc-max", Input).value = "100"
        self.query_one("#score-min", Input).value = ""
        self.query_one("#target-all-check", Checkbox).value = False
        self.notify(
            "Only a few optimal shRNAs matched the default filters "
            "(efficacy score ≥0, GC 30–60%, all targets) — filters have "
            "been relaxed so more candidates are shown.",
            title="Few optimal shRNAs found",
            severity="warning",
            timeout=10,
        )

    def _filtered(self) -> list[Candidate]:
        """Candidates passing the live filter bar, in the current sort order.

        The predicates come from `engine.rejection_reason` — the same function
        `score_one` uses — so this view and a `knockknock run` with equivalent
        settings can never disagree about which candidates pass.
        """
        params = self._export_params()
        results = [
            c
            for c in self._all_candidates
            if rejection_reason(
                c.gc_pct,
                c.efficacy_score,
                c.antisense,
                params,
                targets=c.targets,
            )
            is None
        ]

        results.sort(
            key=lambda c: candidate_sort_key(c, self._sort_by_dg), reverse=True
        )
        return results

    def _read_filter_params(self) -> AnalysisParams:
        """Read the live fields through the engine's shared validator."""
        length = (
            len(self._all_candidates[0].antisense)
            if self._all_candidates
            else AnalysisParams().length
        )
        gc_min_text = self.query_one("#gc-min", Input).value.strip()
        gc_max_text = self.query_one("#gc-max", Input).value.strip()
        score_text = self.query_one("#score-min", Input).value.strip()
        try:
            return validate_analysis_params(
                length=length,
                gc_min=float(gc_min_text) if gc_min_text else 0.0,
                gc_max=float(gc_max_text) if gc_max_text else 100.0,
                min_efficacy_score=float(score_text) if score_text else None,
                sort_by_dg=self._sort_by_dg,
                remove_motifs=self._motif_filter,
                multitarget=self._total_mrnas if self._target_all_only else 1,
                rna_output=self._rna_output,
            )
        except ValueError as exc:
            raise AnalysisParameterError(str(exc)) from exc

    def _refresh(self) -> None:
        # `on_input_changed` schedules this on a 250 ms debounce timer, so it can
        # land after the user has pushed another screen (type a GC value, press
        # `o` straight after). `query_one` on the App searches whichever screen
        # is on top, so the widgets below would raise NoMatches and take the app
        # down. Nothing needs refreshing while the table isn't showing — the
        # explore screen re-reads state from `self` when it comes back.
        try:
            table = self.query_one("#results-table", DataTable)
        except NoMatches:
            return

        candidates = self._filtered()
        table.clear()

        for rank, c in enumerate(candidates, 1):
            mark = "✓" if c.name in self._selected else " "
            antisense, efficacy_score, delta_g_score, gc_pct, dg, targets = (
                format_candidate_cells(c, self._total_mrnas, self._rna_output)
            )
            table.add_row(
                mark,
                str(rank),
                antisense,
                efficacy_score,
                delta_g_score,
                gc_pct,
                dg,
                targets,
                str(c.start_pos),
                c.mrna_id,
                key=c.name,
            )

        n_total = len(self._all_candidates)
        n_shown = len(candidates)
        sel_n = sum(1 for c in candidates if c.name in self._selected)

        settings = []
        if self._motif_filter:
            settings.append("filter low complexity motifs")
        if self._target_all_only:
            settings.append("target all")
        settings.append(f"sort {'ΔG' if self._sort_by_dg else 'efficacy'}")
        filters_text = (
            f"GC {self._gc_min:.0f}–{self._gc_max:.0f}%  |  " + ", ".join(settings)
        )
        if self._filter_error:
            filters_text += f"  |  [red]Filter error: {self._filter_error}[/red]"

        self.query_one("#summary-input", Static).update(
            f"Input: {self._input_path}  |  mRNAs: {self._total_mrnas}"
        )
        self.query_one("#summary-filters", Static).update(filters_text)
        self.query_one("#summary-counts", Static).update(
            f"Showing {n_shown:,} / {n_total:,}  |  Selected: {sel_n:,}"
        )

    def action_toggle_sort(self) -> None:
        self._sort_by_dg = not self._sort_by_dg
        self.query_one("#sort-check", Checkbox).value = self._sort_by_dg
        self._refresh()

    def action_toggle_motifs(self) -> None:
        self._motif_filter = not self._motif_filter
        self.query_one("#motif-check", Checkbox).value = self._motif_filter
        self._refresh()

    def action_toggle_target_all(self) -> None:
        self._target_all_only = not self._target_all_only
        self.query_one("#target-all-check", Checkbox).value = self._target_all_only
        self._refresh()

    def action_focus_gc_min(self) -> None:
        self.query_one("#gc-min", Input).focus()

    def action_focus_gc_max(self) -> None:
        self.query_one("#gc-max", Input).focus()

    def action_focus_score(self) -> None:
        self.query_one("#score-min", Input).focus()

    def action_unfocus(self) -> None:
        """Escape out of a filter-bar field back to the table, without quitting.

        Lets `space`/`o`/etc. work again immediately after typing in a
        filter, instead of having to click the table or tab away first.
        """
        if isinstance(self.focused, (Input, Checkbox)):
            self.query_one("#results-table", DataTable).focus()

    def action_toggle_selection(self) -> None:
        """Toggle the checkmark on the currently-highlighted row."""
        table = self.query_one("#results-table", DataTable)
        if table.row_count == 0:
            return
        try:
            row_key = table.coordinate_to_cell_key(table.cursor_coordinate).row_key
        except CellDoesNotExist:
            # Cursor is parked where a row used to be (filters just narrowed
            # the view). Nothing to toggle.
            return
        name = row_key.value
        if name is None:
            return
        self._toggle_selected_name(name)

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        """Toggle the checkmark on double-click (or Enter) of a row."""
        name = event.row_key.value
        if name is None:
            return
        self._toggle_selected_name(name)

    def _toggle_selected_name(self, name: str) -> None:
        if name in self._selected:
            self._selected.remove(name)
        else:
            self._selected.add(name)
        table = self.query_one("#results-table", DataTable)
        cursor_row = table.cursor_coordinate.row
        self._refresh()
        table = self.query_one("#results-table", DataTable)
        if cursor_row < table.row_count:
            table.move_cursor(row=cursor_row)

    def action_select_all_filtered(self) -> None:
        # `_selected` is a set, so re-adding names is a no-op: 'a' can never
        # double-count a candidate in the "Selected: N" summary.
        for c in self._filtered():
            self._selected.add(c.name)
        self._refresh()

    def action_clear_selection(self) -> None:
        self._selected.clear()
        self._refresh()

    def action_select_top_n(self) -> None:
        table = self.query_one("#results-table", DataTable)
        if table.row_count == 0:
            self.notify("No candidates to select from.", severity="warning")
            return
        self.push_screen(TopNModal(), callback=self._on_top_n_selected)

    def _on_top_n_selected(self, n: int | None) -> None:
        if n is None:
            return
        # Top-N *replaces* the current selection: clear any prior picks first,
        # then check the top-ranked rows of the current view.
        self._selected = {c.name for c in self._filtered()[:n]}
        self._refresh()

    def action_help(self) -> None:
        self.push_screen(HelpScreen())

    def action_open_file(self) -> None:
        """Hand control back to the file picker so another input can be loaded.

        Exiting is the mechanism, not the intent: `launch` runs the picker and
        the explorer as two sibling `App`s in one loop, so returning to the
        picker means ending this app with `open_requested` set. Any pushed
        screen (oligo designer, export chooser) goes with it.
        """
        self.open_requested = True
        self.exit()

    def action_dismiss_banner(self) -> None:
        try:
            banner = self.query_one("#onboarding-banner", Static)
        except NoMatches:
            return  # already dismissed, or never shown
        banner.remove()
        state = _load_state()
        state["seen_banner"] = True
        _save_state(state)

    def _selected_or_filtered_candidates(self) -> list[Candidate]:
        """Selected candidates if any are checked, else everything in the
        current filtered view."""
        filtered = self._filtered()
        if not self._selected:
            return filtered
        candidates = [c for c in filtered if c.name in self._selected]
        if not candidates:
            # The user has selections but none are currently visible
            # (e.g. they tightened filters). Fall back to *all* selected.
            names_by_name = {c.name: c for c in self._all_candidates}
            candidates = [names_by_name[n] for n in self._selected if n in names_by_name]
        return candidates

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "design-oligos-btn":
            self.action_design_oligos()

    def action_design_oligos(self) -> None:
        candidates = self._selected_or_filtered_candidates()
        if not candidates:
            self.notify("No candidates to design oligos for.", severity="warning")
            return
        self.push_screen(
            OligoScreen(
                candidates,
                input_path=self._input_path,
                params=self._export_params(),
                total_mrnas=self._total_mrnas,
            )
        )

    def action_blast_selected(self) -> None:
        """Open NCBI Nucleotide BLAST in the browser, prefilled with the
        selected guide sequences (or the current filtered view if nothing
        is checked)."""
        candidates = self._selected_or_filtered_candidates()
        if not candidates:
            self.notify("No candidates to BLAST.", severity="warning")
            return

        truncated = len(candidates) > _BLAST_MAX_SEQS
        candidates = candidates[:_BLAST_MAX_SEQS]
        fasta = "\n".join(f">{c.name}\n{c.antisense}" for c in candidates)
        url = (
            f"{_BLAST_URL}?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC="
            f"&LINK_LOC=blasthome&BLAST_PROGRAMS=blastn"
            f"&QUERY={quote(fasta, safe='')}"
        )
        webbrowser.open(url)

        message = f"Opened NCBI BLAST with {len(candidates)} guide sequence(s)."
        if truncated:
            message += (
                f" Capped at {_BLAST_MAX_SEQS} — narrow your selection to BLAST more."
            )
        self.notify(message, timeout=6)

    def _export_params(self) -> AnalysisParams:
        """Parameters describing the *current view*, so the provenance block in
        an export records the filters actually applied in the TUI rather than
        the wide-open ones the explore session was started with."""
        length = (
            len(self._all_candidates[0].antisense)
            if self._all_candidates
            else AnalysisParams().length
        )
        return validate_analysis_params(
            length=length,
            gc_min=self._gc_min,
            gc_max=self._gc_max,
            min_efficacy_score=self._score_min,
            sort_by_dg=self._sort_by_dg,
            remove_motifs=self._motif_filter,
            multitarget=self._total_mrnas if self._target_all_only else 1,
            rna_output=self._rna_output,
        )

    def action_export_view(self) -> None:
        if not self._filtered():
            self.notify("No candidates to export.", severity="warning")
            return
        self.push_screen(
            ExportFormatModal(CANDIDATE_FORMATS, title="Export the current view as:"),
            callback=self._on_export_format_chosen,
        )

    def _on_export_format_chosen(self, fmt: str | None) -> None:
        if fmt is None:
            return

        candidates = self._filtered()
        params = self._export_params()
        ext = {f.value: f.ext for f in CANDIDATE_FORMATS}[fmt]
        out_path = resolve_export_path(
            ext,
            reference_path=self._input_path,
            stem="knockknock_export",
            label="candidates",
        )

        def _write() -> None:
            write_candidates_export(
                fmt,
                candidates,
                params,
                out_path,
                input_path=self._input_path,
                total_mrnas=self._total_mrnas,
            )

        export_and_notify(
            self, out_path, _write, count=len(candidates), noun="candidates"
        )

    def on_input_changed(self, event: Input.Changed) -> None:
        try:
            params = self._read_filter_params()
        except AnalysisParameterError as exc:
            # Keep filtering on the last valid settings while showing exactly
            # why the in-progress text is not accepted.
            self._filter_error = str(exc)
        else:
            self._gc_min = params.gc_min
            self._gc_max = params.gc_max
            self._score_min = params.min_efficacy_score
            self._filter_error = None
        # Debounce: filtering a few-thousand-candidate list on every keystroke
        # makes typing feel laggy. Coalesce successive edits into a single
        # refresh after 250 ms of idle.
        if self._refresh_timer is not None:
            self._refresh_timer.stop()
        self._refresh_timer = self.set_timer(0.25, self._refresh)

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        if event.checkbox.id == "motif-check":
            self._motif_filter = event.checkbox.value
        elif event.checkbox.id == "sort-check":
            self._sort_by_dg = event.checkbox.value
        elif event.checkbox.id == "target-all-check":
            self._target_all_only = event.checkbox.value
        self._refresh()


_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn", ".txt"}


def tree_root_for(path: Path | None) -> Path:
    """The directory a picker tree should be rooted at to show ``path``.

    With no path to show (the picker's first launch), root at the current
    working directory — that's usually closer to the data than home is.
    `DirectoryTree` can't natively step above its own root, so this alone
    would wall the user in; `_FastaDirectoryTree.action_go_up` (bound to
    `backspace`) is what re-roots the tree one level up at a time, all the
    way to the filesystem root, so home and everything beside the cwd stay
    reachable.

    With a path to show (``Ctrl+O`` reopening the picker on a previously
    loaded file), root at home instead and expand down to it. A tree can't
    be navigated above its own root, so rooting it at the folder the file
    happens to sit in would hide everything around it — you couldn't step up
    to a sibling directory, and the root row would read as a bare folder
    name with no clue where it is. Rooting at home and expanding down
    instead leaves the whole path on screen and the rest of the filesystem
    reachable.

    Files outside home (an external volume, ``/tmp``) fall back to their
    filesystem anchor, which is the nearest root that can still contain them.
    """
    if path is None:
        return Path.cwd()
    home = Path.home()
    if path == home or home in path.parents:
        return home
    return Path(path.anchor or home)


class _ParentDirEntry(DirEntry):
    """Marks the synthetic ".." row so it re-roots the tree instead of being
    treated as a real, expandable child of the root."""


class _FastaDirectoryTree(DirectoryTree):
    """DirectoryTree that only lists directories and FASTA-like files.

    ``highlight`` opens the tree with one file already under the cursor and
    every directory above it expanded — how `PickerApp` puts the
    previously-loaded input back where the user left it, in context.

    The root's own row is hidden (``show_root``) so the ".." shortcut
    `_populate_node` adds can sit at the very top of the listing, line 0,
    rather than one line below a redundant "current folder" row. The
    current folder is shown instead by ``on_root_changed``, which
    `PickerApp` uses to keep a breadcrumb `Static` above the tree in sync.
    """

    _FASTA_SUFFIXES = _FASTA_SUFFIXES

    BINDINGS = [
        Binding("backspace", "go_up", "Up a folder"),
        Binding("enter", "select_cursor", "Open", show=True),
    ]

    def __init__(
        self,
        path,
        *,
        highlight: Path | None = None,
        on_root_changed: Callable[[Path], None] | None = None,
        **kwargs,
    ):
        super().__init__(path, **kwargs)
        self.show_root = False
        self._highlight = highlight.resolve() if highlight is not None else None
        #: Where the tree started. The highlight is only restored while the
        #: tree still shows this root, so navigating elsewhere (dropping a
        #: directory onto the window) doesn't drag the cursor back.
        self._initial_root = Path(self.path)
        self._on_root_changed = on_root_changed

    def filter_paths(self, paths):  # type: ignore[override]
        out = []
        for p in paths:
            if p.name.startswith("."):
                continue
            if p.is_dir():
                out.append(p)
            elif p.suffix.lower() in self._FASTA_SUFFIXES:
                out.append(p)
        return out

    def _populate_node(self, node, content):  # type: ignore[override]
        super()._populate_node(node, content)
        # Only the root gets a ".." row: it's the only level `action_go_up`
        # can actually act on, since every other row is just a nested
        # directory that already collapses in place.
        if node is not self.root:
            return
        parent = Path(self.path).parent
        if parent == Path(self.path):
            return  # already at the filesystem root — nothing above to show
        node.add("..", data=_ParentDirEntry(parent), allow_expand=False, before=0)

    def render_label(self, node, base_style, style):  # type: ignore[override]
        if isinstance(node.data, _ParentDirEntry):
            return Text("⬆ .. (up one level)", style=base_style + style)
        return super().render_label(node, base_style, style)

    async def _on_tree_node_selected(self, event: Tree.NodeSelected) -> None:
        if isinstance(event.node.data, _ParentDirEntry):
            event.stop()
            self.action_go_up()
            return
        await super()._on_tree_node_selected(event)

    async def on_mount(self) -> None:
        await self._restore_highlight()

    async def watch_path(self) -> None:
        # The base watcher repopulates the tree and then parks the cursor on
        # line 0, so restoring the highlight has to run after it — doing this
        # from the owning app's `on_mount` instead loses the race.
        await super().watch_path()
        await self._restore_highlight()
        if self._highlight is None:
            self._land_past_parent_row()
        if self._on_root_changed is not None:
            self._on_root_changed(Path(self.path))

    def _land_past_parent_row(self) -> None:
        """On a fresh load with nothing to highlight, skip past the synthetic
        ".." row so the cursor starts on real content, not the shortcut."""
        children = self.root.children
        if len(children) > 1 and isinstance(children[0].data, _ParentDirEntry):
            _ = self.last_line
            self.move_cursor(children[1], animate=False)

    async def _restore_highlight(self) -> None:
        """Expand the chain of directories down to ``highlight`` and select it.

        Walks one level at a time because that is how `DirectoryTree` works: a
        directory's children don't exist until it has been loaded, so each
        ``await`` is what makes the next level searchable.

        Giving up part-way is fine and expected — the file may have been
        filtered out by `filter_paths` (a FASTA with an unlisted suffix, named
        on the command line), or a directory above it may be unreadable. The
        tree is then left expanded as far as it got, which is still useful.
        """
        if self._highlight is None or Path(self.path) != self._initial_root:
            return
        try:
            parts = self._highlight.relative_to(self._initial_root).parts
        except ValueError:
            return  # highlight isn't under this root after all

        node: TreeNode[DirEntry] = self.root
        for part in parts:
            await self.reload_node(node)
            node.expand()
            next_node = next(
                (
                    child
                    for child in node.children
                    if child.data is not None and child.data.path.name == part
                ),
                None,
            )
            if next_node is None:
                return
            node = next_node

        # The loop above only expands nodes on the way *to* the highlight,
        # not the highlight itself. That's correct when it's a file (nothing
        # to expand), but when `action_go_up` feeds in the folder just left,
        # leaving it collapsed would hide the very listing the user was
        # just looking at.
        if node.allow_expand:
            await self.reload_node(node)
            node.expand()

        # `move_cursor` positions by the node's cached line number, and every
        # `expand` above invalidated that cache. Reading `last_line` is what
        # rebuilds it; without this the cursor is set to -1 and snaps back to
        # the root. `move_cursor` scrolls the node into view itself.
        _ = self.last_line
        self.move_cursor(node, animate=False)

    def action_go_up(self) -> None:
        """Re-root the tree one level up, with the folder we came from selected.

        `DirectoryTree` has no notion of a parent — it can only show its own
        root and below — so "going up" means discarding the current root and
        starting over one level higher. Feeding the old root back in as
        ``_highlight`` reuses `_restore_highlight` to land the cursor right
        back where it was, instead of dropping the user at the top of an
        unfamiliar listing.
        """
        current_root = Path(self.path)
        parent = current_root.parent
        if parent == current_root:
            return  # already at the filesystem root
        self._highlight = current_root
        self._initial_root = parent
        self.path = str(parent)


def _parse_dropped_path(text: str) -> Path | None:
    """Best-effort extraction of a file path from dropped/pasted terminal text.

    Dragging a file onto Terminal.app/iTerm2 inserts its path as pasted text,
    either backslash-escaped or wrapped in quotes, so a plain ``Path(text)``
    often fails to resolve. ``shlex.split`` undoes that escaping.
    """
    text = text.strip()
    if not text:
        return None
    try:
        parts = shlex.split(text)
    except ValueError:
        parts = [text]
    if not parts:
        return None
    path = Path(parts[0]).expanduser()
    return path if path.exists() else None


class ProcessingScreen(Screen):
    """Shown while the chosen FASTA file is parsed and shRNAs are designed.

    Parsing and scoring can take several seconds on larger inputs; without
    this screen the terminal briefly looks like it dropped back to the shell
    (or crashed) between file selection and the results table appearing.
    """

    DEFAULT_CSS = """
    ProcessingScreen {
        align: center middle;
    }
    ProcessingScreen > Vertical {
        width: 50;
        height: auto;
        align: center middle;
    }
    ProcessingScreen LoadingIndicator {
        width: 100%;
        height: 3;
    }
    ProcessingScreen Static {
        width: 100%;
        margin-top: 1;
        text-align: center;
    }
    """

    def __init__(self, message: str):
        super().__init__()
        self._message = message

    def compose(self) -> ComposeResult:
        with Vertical():
            yield LoadingIndicator()
            yield Static(self._message)


class PickerApp(App):
    """Standalone file-picker app used when ``knockknock`` is invoked without args.

    Offers two ways in: pick (or drag in) a FASTA file, or paste a raw
    sequence/FASTA into the text box below the tree. Either way, a
    :class:`ProcessingScreen` is shown while the input is parsed and shRNA
    candidates are designed in a background worker, so the caller can hand
    the finished results straight to :class:`ExploreApp`. Results are
    stashed on ``self.candidates``/``self.total_mrnas``/``self.input_path``
    after ``run()`` returns; ``self.candidates`` is still ``None`` if the user
    quit without choosing anything.

    **Nothing here quits the app on its own.** An unreadable file, a broken
    FASTA, or an input that yields no candidates is reported as a toast and
    leaves the picker up, because the useful next step in every one of those
    cases is picking a different file — and this is also where `ExploreApp`'s
    ``ctrl+o`` lands, so quitting would throw away a loaded table.

    ``initial_path``/``initial_text`` pre-seed the picker with whatever was
    open before: the tree opens in that file's directory with it under the
    cursor, or the paste box comes back with the pasted sequence in it.
    """

    TITLE = "KnockKnock"
    SUB_TITLE = "choose a FASTA file, or paste a sequence"

    CSS = """
    Screen {
        layout: vertical;
    }
    #picker-help {
        height: auto;
        padding: 1 2;
        background: $surface;
        border-bottom: solid $primary;
    }
    #picker-path {
        height: auto;
        padding: 0 2;
        background: $surface;
        color: $text-muted;
    }
    DirectoryTree {
        height: 1fr;
    }
    #paste-section {
        height: auto;
        padding: 1 2;
        background: $boost;
        border-top: solid $primary;
    }
    #paste-section > Static {
        margin-bottom: 1;
    }
    #paste-area {
        height: 6;
        border: solid $primary-background;
    }
    #paste-submit {
        margin-top: 1;
        width: auto;
        min-width: 30;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("escape", "unfocus", "Unfocus", show=False),
    ]

    def __init__(
        self,
        initial_path: Path | None = None,
        initial_text: str | None = None,
    ):
        super().__init__()
        self._initial_path = initial_path
        self._initial_text = initial_text
        #: What the user picked, in the two shapes `launch` needs to seed the
        #: picker again next time. Exactly one is ever set: a file open clears
        #: the pasted text and vice versa, so the seed is never ambiguous.
        self.selected_path: Path | None = None
        self.source_text: str | None = None
        self.candidates: list[Candidate] | None = None
        self.total_mrnas: int = 1
        self.input_path: str = ""
        self.warnings: list[str] = []

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        yield Static(
            "[b]KnockKnock[/b]  —  pick a FASTA file, or drag one in "
            "(Enter to open, Backspace for parent folder, q to quit)",
            id="picker-help",
        )
        yield Static(str(tree_root_for(self._initial_path)), id="picker-path")
        yield _FastaDirectoryTree(
            str(tree_root_for(self._initial_path)),
            highlight=self._initial_path,
            on_root_changed=self._set_path_label,
            id="picker-tree",
        )
        with Vertical(id="paste-section"):
            yield Static(
                "[dim]— or paste a sequence / FASTA below, then click Design —[/dim]"
            )
            yield TextArea(id="paste-area")
            yield Button("Design from Pasted Text", id="paste-submit", variant="primary")
        yield Footer()

    def check_action(self, action: str, parameters: tuple) -> bool | None:
        # Don't let 'q' quit out from under someone mid-paste — FASTA
        # headers routinely contain 'q' (e.g. UniProt accessions like Q9Y6K9).
        if action == "quit" and isinstance(self.focused, TextArea):
            return False
        return True

    def get_key_display(self, binding: Binding) -> str:
        return _ctrl_key_display(binding) or super().get_key_display(binding)

    def action_unfocus(self) -> None:
        """Escape out of the paste box back to the tree, without quitting."""
        if isinstance(self.focused, TextArea):
            self.query_one("#picker-tree", _FastaDirectoryTree).focus()

    def on_mount(self) -> None:
        """Put the previously pasted sequence back in the box, if there was one.

        The file-tree half of the seeding is `_FastaDirectoryTree`'s own job —
        it has to outlast the cursor reset its base class does on load.
        """
        if self._initial_text is not None:
            area = self.query_one("#paste-area", TextArea)
            area.text = self._initial_text
            area.focus()

    def _set_path_label(self, path: Path) -> None:
        """Keep the breadcrumb above the tree in sync with its (hidden) root.

        The tree's own root row is hidden so the ".." shortcut can sit at
        line 0 instead of one line below a redundant "current folder" row —
        this label is what shows that folder instead.
        """
        self.query_one("#picker-path", Static).update(str(path))

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        self._open_file(Path(event.path).resolve())

    def on_paste(self, event: events.Paste) -> None:
        """Handle a file dragged onto the terminal window.

        Terminals deliver a drag-and-drop as a paste of the file's path
        rather than a distinct drop event, so we sniff pasted text for a
        path and open it directly if it looks like one.
        """
        path = _parse_dropped_path(event.text)
        if path is None:
            return
        if path.is_dir():
            self.query_one("#picker-tree", _FastaDirectoryTree).path = str(path)
        elif path.is_file():
            self._open_file(path.resolve())

    def _open_file(self, path: Path) -> None:
        self.selected_path = path
        self.source_text = None
        self.push_screen(ProcessingScreen(f"Parsing {path.name} and designing shRNAs…"))
        # Reading happens inside the worker, not here: a file picked from the
        # tree can turn out to be binary or non-UTF-8 (".txt" is an accepted
        # suffix), and that decode error must surface as a message rather than
        # taking the app down mid-event-handler.
        self._process_input(path=path, source_label=str(path), error_noun=str(path))

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "paste-submit":
            self._submit_pasted_text()

    def _submit_pasted_text(self) -> None:
        text = self.query_one("#paste-area", TextArea).text
        if not text.strip():
            self.notify("Paste a sequence or FASTA first.", severity="warning")
            return
        self.selected_path = None
        self.source_text = text
        self.push_screen(ProcessingScreen("Parsing pasted input and designing shRNAs…"))
        self._process_input(text=text, source_label="-", error_noun="pasted input")

    @work(thread=True)
    def _process_input(
        self,
        *,
        source_label: str,
        error_noun: str,
        text: str | None = None,
        path: Path | None = None,
    ) -> None:
        """Load ``text`` (or the contents of ``path``) and design shRNAs for it.

        Shared by the file picker and the paste box, and — via
        `seqinput.load_sequences` — by the CLI too, so alphabet validation,
        duplicate-header rejection and U→T normalisation are identical however
        the sequence arrived.
        """
        try:
            if path is not None:
                text = read_sequence_file(path)
            mrnas, warnings = load_sequences(text or "")
            params = AnalysisParams()
            candidates_by_mrna = make_shrnas(mrnas, params)
            candidates, stats = process_shrnas(candidates_by_mrna, mrnas, params)
        except SequenceInputError as exc:
            self.call_from_thread(self._processing_failed, str(exc))
            return
        except Exception as exc:
            self.call_from_thread(
                self._processing_failed, f"Could not process {error_noun}: {exc}"
            )
            return

        self.warnings = warnings
        # The picker runs the engine with default (wide-open) parameters, so an
        # empty result is not something the explorer's filters could recover
        # from — report it here instead of handing over an empty table.
        if not candidates:
            self.call_from_thread(
                self._processing_failed,
                no_candidates_message(stats, params),
                "warning",
            )
            return

        self.candidates = candidates
        self.total_mrnas = len(mrnas)
        self.input_path = source_label
        self.call_from_thread(self.exit)

    def _processing_failed(
        self,
        message: str,
        severity: Literal["information", "warning", "error"] = "error",
    ) -> None:
        """Report a bad input and hand the picker back, rather than quitting.

        ``severity`` separates "this input is broken" from "this input is fine
        but has no candidates in it"; both leave the user where they can pick
        something else.
        """
        if isinstance(self.screen, ProcessingScreen):
            self.pop_screen()
        self.notify(message, severity=severity, timeout=12)


def launch(
    candidates: list[Candidate] | None = None,
    total_mrnas: int = 1,
    input_path: str = "",
    rna_output: bool = False,
    warnings: list[str] | None = None,
) -> None:
    """Launch the interactive explore TUI.

    If ``candidates`` is ``None`` we open a file picker first, parse the chosen
    FASTA, run the engine with default parameters, and then hand off to the
    explore screen. Otherwise the candidates are shown directly.

    ``warnings`` are input notes to surface once the explorer is up. They are
    handed to `ExploreApp` rather than printed here, because the app clears
    the screen the moment it starts.

    The loop is what makes `ExploreApp`'s ``ctrl+o`` work. Textual apps aren't
    nestable — an `App` can push screens but not another `App` — so "open a
    different file" can't be a screen swap inside the explorer. Instead the
    explorer ends with `ExploreApp.open_requested` set, and we come back round
    to the picker, seeded with whatever was open so the user starts next to it
    rather than back at their home directory.
    """
    warnings = warnings or []
    seed_path: Path | None = None
    seed_text: str | None = None
    if candidates is not None and input_path and input_path != "-":
        # Launched straight into the explorer (`run --tui`): the only clue
        # about what to seed the picker with is the input path we were given.
        given = Path(input_path).expanduser()
        if given.is_file():
            seed_path = given.resolve()

    while True:
        if candidates is None:
            picker = PickerApp(initial_path=seed_path, initial_text=seed_text)
            picker.run()
            if picker.candidates is None:
                return  # user quit without selecting anything
            candidates = picker.candidates
            total_mrnas = picker.total_mrnas
            input_path = picker.input_path
            warnings = picker.warnings
            seed_path = picker.selected_path
            seed_text = picker.source_text

        app = ExploreApp(
            candidates,
            total_mrnas,
            input_path,
            rna_output=rna_output,
            warnings=warnings,
        )
        app.run()
        if not app.open_requested:
            return
        candidates = None
