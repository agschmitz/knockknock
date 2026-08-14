"""Shared base for small prompt modals (`TopNModal`, `SaveBackboneModal`,
`ExportFormatModal`) — just the part of their behavior that's identical
across all three: Escape always cancels by dismissing with `None`."""

from typing import TypeVar

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal
from textual.screen import ModalScreen
from textual.widgets import Button, Static

ResultT = TypeVar("ResultT")


class PromptModal(ModalScreen[ResultT]):
    """A modal screen that dismisses with `None` when the user presses Escape.

    Subclasses provide their own `compose()`/`DEFAULT_CSS`/other bindings;
    inheriting this keeps the escape-cancels contract from drifting between
    dialogs that each used to redeclare it.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    def action_cancel(self) -> None:
        self.dismiss(None)


class ConfirmModal(PromptModal[bool]):
    """Yes/no confirmation for a destructive action.

    Dismisses with ``True`` only on an explicit confirmation; Escape and the
    Cancel button both give ``False`` via `PromptModal`'s contract, so a
    stray keypress can never delete anything.
    """

    DEFAULT_CSS = """
    ConfirmModal {
        align: center middle;
    }
    ConfirmModal > Container {
        width: auto;
        min-width: 50;
        max-width: 70;
        height: auto;
        background: $surface;
        border: thick $error;
        padding: 1 2;
    }
    ConfirmModal Horizontal {
        height: auto;
        margin-top: 1;
        align: center middle;
    }
    ConfirmModal Button {
        margin: 0 1;
    }
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("n", "cancel", "No"),
        Binding("y", "confirm", "Yes"),
    ]

    def __init__(self, question: str, *, detail: str = "") -> None:
        super().__init__()
        self._question = question
        self._detail = detail

    def compose(self) -> ComposeResult:
        with Container():
            yield Static(self._question)
            if self._detail:
                yield Static(f"[dim]{self._detail}[/dim]")
            with Horizontal():
                yield Button("Delete", id="confirm-yes", variant="error")
                yield Button("Cancel", id="confirm-no", variant="primary")
            yield Static("[dim]y to confirm, n or Esc to cancel[/dim]")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.dismiss(event.button.id == "confirm-yes")

    def action_confirm(self) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)
