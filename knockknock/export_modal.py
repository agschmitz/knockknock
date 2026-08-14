"""Modal export-format chooser shared by the explore and oligo TUI screens.

The CLI has always been able to write several formats (``--format``); the TUI
used to hard-code one per screen (CSV for the candidate table, TXT for the
oligo designer). This modal gives the interactive UI the same choice, and
lives in its own module so both screens can import it without ``oligo_tui``
having to import from ``tui`` (which imports ``oligo_tui``).
"""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from textual.app import ComposeResult
from textual.containers import Container
from textual.widgets import OptionList, Static
from textual.widgets.option_list import Option

from knockknock.modals import PromptModal


@dataclass(frozen=True)
class ExportFormat:
    """One row in the chooser.

    ``value`` is the identifier handed back to the calling screen; it matches
    the CLI's ``--format`` spelling so the two stay in step.
    """

    value: str
    label: str
    ext: str
    hint: str


#: Formats offered for the candidate table (explore screen). Same set — and
#: same spellings — as ``knockknock --format``, minus ``table``, which is what
#: the screen already shows.
CANDIDATE_FORMATS: tuple[ExportFormat, ...] = (
    ExportFormat("csv", "CSV", ".csv", "one row per guide, with run metadata"),
    ExportFormat(
        "xlsx", "Excel", ".xlsx", "coloured efficacy/ΔG score scales + Run Info sheet"
    ),
    ExportFormat("json", "JSON", ".json", "structured, for scripts and pipelines"),
    ExportFormat("fasta", "FASTA", ".fasta", "guide sequences for BLAST etc."),
)

#: Formats offered for assembled oligo pairs (oligo designer screen).
OLIGO_FORMATS: tuple[ExportFormat, ...] = (
    ExportFormat("fasta", "FASTA", ".fasta", "forward/reverse oligo pairs"),
    ExportFormat("csv", "CSV", ".csv", "one row per pair, for ordering sheets"),
    ExportFormat("json", "JSON", ".json", "structured, for scripts and pipelines"),
)


def export_and_notify(
    screen, out_path: Path, write: Callable[[], None], *, count: int, noun: str
) -> None:
    """Run a format-specific ``write()``, then notify success/failure.

    `tui.py`'s and `oligo_tui.py`'s export-modal callbacks both resolve a
    path, write the chosen format, and notify the same way — this is that
    shared write/notify tail so the two can't drift on the message wording
    or on catching the same `OSError`.
    """
    try:
        write()
    except OSError as exc:
        screen.notify(f"Could not write {out_path}: {exc}", severity="error")
        return
    screen.notify(f"Exported {count} {noun} to [bold]{out_path}[/bold]", timeout=5)


class ExportFormatModal(PromptModal[str | None]):
    """Ask which format to export as. Dismisses with an ``ExportFormat.value``,
    or ``None`` if the user cancelled."""

    DEFAULT_CSS = """
    ExportFormatModal {
        align: center middle;
    }
    ExportFormatModal > Container {
        width: auto;
        min-width: 56;
        height: auto;
        background: $surface;
        border: thick $primary;
        padding: 1 2;
    }
    ExportFormatModal OptionList {
        width: auto;
        height: auto;
        max-height: 12;
        border: none;
        background: $surface;
        margin: 1 0;
    }
    """

    def __init__(
        self,
        formats: tuple[ExportFormat, ...] = CANDIDATE_FORMATS,
        *,
        title: str = "Export as:",
    ) -> None:
        super().__init__()
        self._formats = tuple(formats)
        self._title = title

    def compose(self) -> ComposeResult:
        # Pad the label column so the extension/hint columns line up whatever
        # mix of formats a screen offers.
        width = max((len(f.label) for f in self._formats), default=0)
        with Container():
            yield Static(self._title)
            yield OptionList(
                *(
                    Option(
                        f"[b]{i}[/b]  {fmt.label:<{width}}  "
                        f"[dim]{fmt.ext} — {fmt.hint}[/dim]",
                        id=fmt.value,
                    )
                    for i, fmt in enumerate(self._formats, 1)
                ),
                id="export-format-list",
            )
            yield Static("[dim]1–9 to pick directly, ↑↓ then Enter, Esc to cancel[/dim]")

    def on_mount(self) -> None:
        self.query_one("#export-format-list", OptionList).focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.dismiss(event.option.id)

    def on_key(self, event) -> None:
        """Number keys jump straight to a format, so a two-keystroke export
        (``e`` then ``2``) stays as fast as the old single-format binding."""
        char = event.character
        if char is None or not char.isdigit():
            return
        index = int(char) - 1
        if 0 <= index < len(self._formats):
            event.stop()
            event.prevent_default()
            self.dismiss(self._formats[index].value)
