"""`TwoRowFooter` — a `Footer` that wraps bindings onto two lines.

Textual's `Footer` docks a single `height: 1` row and silently clips
whatever doesn't fit (mouse-wheel-scrollable, but with no visible
indication there's more). The explorer and oligo designer screens each
bind well over a dozen keys, so a single row runs out of width on any
normal terminal. Splitting them across two rows instead of shrinking or
hiding any of them is the fix — but *where* to split matters: bindings
are declared in each screen's `BINDINGS` in related clusters (filters,
selection, actions, ...), and `row_break` names the action that should
start the second row, so that clustering survives onto the two lines
instead of getting cut at whatever the midpoint happens to be. The
command-palette key (when enabled) always lands at the end of row two.

Grouped bindings (`Binding(group=...)`) aren't used anywhere in this app,
so that part of `Footer.compose` isn't reimplemented.
"""

from __future__ import annotations

from math import ceil

from textual.app import ComposeResult
from textual.containers import HorizontalGroup
from textual.widgets import Footer
from textual.widgets._footer import FooterKey


class TwoRowFooter(Footer):
    """Footer variant that lays active bindings out across two rows.

    `row_break` is the action name of the first binding that belongs on
    the second row; every binding before it (in `Screen.active_bindings`
    order, which follows `BINDINGS` declaration order) goes on the first.
    """

    DEFAULT_CSS = """
    TwoRowFooter {
        height: 2;
        layout: vertical;
    }
    """

    def __init__(self, row_break: str, **kwargs) -> None:
        super().__init__(**kwargs)
        self._row_break = row_break

    def compose(self) -> ComposeResult:
        if not self._bindings_ready:
            return
        active_bindings = self.screen.active_bindings
        bindings = [
            (binding, enabled, tooltip)
            for (_, binding, enabled, tooltip) in active_bindings.values()
            if binding.show
        ]
        if self.show_command_palette and self.app.ENABLE_COMMAND_PALETTE:
            try:
                bindings.append(active_bindings[self.app.COMMAND_PALETTE_BINDING][1:])
            except KeyError:
                pass

        split = next(
            (i for i, (binding, _, _) in enumerate(bindings) if binding.action == self._row_break),
            ceil(len(bindings) / 2),
        )
        for row_bindings in (bindings[:split], bindings[split:]):
            with HorizontalGroup():
                for binding, enabled, tooltip in row_bindings:
                    yield FooterKey(
                        binding.key,
                        self.app.get_key_display(binding),
                        binding.description,
                        binding.action,
                        disabled=not enabled,
                        tooltip=tooltip,
                    ).data_bind(compact=Footer.compact)
