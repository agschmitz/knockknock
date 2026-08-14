"""Modal help overlay used by both the explore and oligo TUI screens."""

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Container, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import Footer, Static

_HELP_MARKDOWN = """\
[b]KnockKnock — help[/b]

[b]Columns[/b]
  [cyan]Efficacy Score[/cyan]  Sum of Matveeva position-specific coefficients.
                Higher is better; a rough rule of thumb is ≥ 0.3.
  [cyan]ΔG Score[/cyan]    Binned ΔG penalty (Xia 1998 nearest-neighbor).
                Higher = more favourable strand selection.
  [cyan]GC%[/cyan]         GC content of the guide. Typical good range: 30–60 %.
  [cyan]ΔG[/cyan]          ΔG of the guide-target RNA duplex (kcal/mol).
                Typical good range: -33 ± 7. More negative = more stable.
  [cyan]Targets[/cyan]     How many of the input mRNAs share this guide.

[b]Explore screen[/b]
  space         Toggle selection on the current row
  a             Select all currently-filtered candidates
  Shift+A       Clear all selections
  n             Select top N candidates from the current view
  s             Toggle sort (efficacy ↔ ΔG)
  m             Toggle low complexity motif filter
  t             Toggle "target all inputs" filter
  Shift+T       Focus the efficacy score input
  g / Shift+G   Focus the GC-min / GC-max input
  o             Design oligos (uses selection, or all-filtered if none)
  b             BLAST selection at NCBI (opens browser, capped at 25 seqs)
  e             Export the current view (choose CSV / Excel / JSON / FASTA)
  Ctrl+O        Open another input (back to the file picker)
  ?             This help
  q             Quit

[b]Oligo screen[/b]
  ← →           Previous / next candidate
  c             Copy the example oligo pair to clipboard
  e             Export all oligo pairs (choose text / CSV / JSON)
  s             Save the current junctions/flanks/loop as a new backbone
  d             Delete the selected backbone (bring bundled ones back with
                --restore-backbones)
  b / Esc       Back to the explore screen
  ?             This help

[b]Export format chooser[/b]
  1–9           Pick a format directly
  ↑↓ + Enter    Pick from the list
  Esc           Cancel without writing a file

[b]Quality presets (CLI)[/b]
  --quality lax        no filters
  --quality standard   GC 20–80, motif filter on
  --quality strict     GC 30–60, motif filter on, min efficacy score ≥ 0

[dim]Press Esc or ? to close.[/dim]
"""


class HelpScreen(ModalScreen):
    """Modal help / cheat-sheet overlay. Press Esc or ? to dismiss."""

    BINDINGS = [
        Binding("escape", "dismiss", "Close"),
        Binding("question_mark", "dismiss", "Close"),
        Binding("q", "dismiss", "Close"),
    ]

    DEFAULT_CSS = """
    HelpScreen {
        align: center middle;
    }
    HelpScreen > Container {
        width: 80%;
        max-width: 90;
        height: auto;
        max-height: 90%;
        background: $surface;
        border: thick $primary;
        padding: 1 2;
    }
    HelpScreen VerticalScroll {
        height: auto;
        max-height: 30;
    }
    HelpScreen Static {
        height: auto;
    }
    """

    def compose(self) -> ComposeResult:
        with Container():
            with VerticalScroll():
                yield Static(_HELP_MARKDOWN)
        yield Footer()
