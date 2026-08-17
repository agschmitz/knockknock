"""TUI tests using Textual's Pilot harness.

These tests don't try to assert on the rendered pixels; they verify behaviour:
state changes after specific key/button presses, mounted widget queries,
and that screens can be opened and closed without errors.
"""

from pathlib import Path

import pytest
from Bio.Seq import Seq
from textual.worker import WorkerCancelled

from knockknock.engine import (
    AnalysisParams,
    make_shrnas,
    process_shrnas,
)
from knockknock.export_modal import CANDIDATE_FORMATS, ExportFormatModal
from knockknock.oligo_tui import OligoScreen
from knockknock.tui import ExploreApp, PickerApp, tree_root_for

pytestmark = pytest.mark.asyncio


@pytest.fixture
def candidates():
    mrnas = [("test", Seq("ATGCGTAGCTAGCTAGCTAGCTGATCGTAGCTAGCTAGCT"))]
    params = AnalysisParams()
    cmap = make_shrnas(mrnas, params)
    results, _ = process_shrnas(cmap, mrnas, params)
    assert results, "test fixture should produce candidates"
    return results


async def test_explore_screen_mounts_and_shows_rows(candidates):
    """The results table is populated on mount, with one row per candidate
    that survives the (wide-open, default) live filters."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        table = app.query_one("#results-table")
        assert table.row_count > 0
        assert table.row_count == len(app._filtered())


async def test_invalid_filter_range_is_rejected_without_changing_live_filter(candidates):
    """The TUI keeps the last valid values and exposes the shared error."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#gc-min").value = "80"
        app.query_one("#gc-max").value = "20"
        await pilot.pause()

        assert app._gc_min == 30.0
        assert app._gc_max == 60.0
        assert app._filter_error is not None
        assert "GC limits" in app._filter_error


async def test_space_toggles_selection(candidates):
    """Space checks/unchecks the row under the cursor in `_selected`."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        assert len(app._selected) == 0
        # Move focus to the data table, then press space.
        table = app.query_one("#results-table")
        table.focus()
        await pilot.press("space")
        await pilot.pause()
        assert len(app._selected) == 1
        # Pressing again toggles it off.
        await pilot.press("space")
        await pilot.pause()
        assert len(app._selected) == 0


async def test_double_click_toggles_selection(candidates):
    """Double-clicking a row checks/unchecks it in `_selected`, same as space."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    # A taller terminal than the default so the results table (below the
    # filter bar and summary panel) actually has rows to click on.
    async with app.run_test(size=(120, 50)) as pilot:
        await pilot.pause()
        assert len(app._selected) == 0
        # Row 1 (below the header at local y=0 and row 0 at local y=2).
        row1_offset = (2, 3)
        await pilot.click("#results-table", offset=row1_offset, times=2)
        await pilot.pause()
        assert len(app._selected) == 1
        # Move the cursor away so the next double-click's first click has to
        # land the cursor on row 1 again before its second click can toggle it
        # (mirrors DataTable's "click to move, click again to select" model).
        await pilot.press("up")
        await pilot.pause()
        await pilot.click("#results-table", offset=row1_offset, times=2)
        await pilot.pause()
        assert len(app._selected) == 0


@pytest.fixture
def multi_mrna_candidates():
    # Both sequences share a common region so some guides will target both,
    # while the unique tail of each sequence produces mRNA-only guides.
    shared = "ATGCGTAGCTAGCTAGCTAGCTGATCGTAGCTAGCTAGCT"
    mrnas = [
        ("mrna_a", Seq(shared + "AAAAACCCCCGGGGGTTTTT")),
        ("mrna_b", Seq(shared + "GGGGGAAAAATTTTTCCCCC")),
    ]
    params = AnalysisParams()
    cmap = make_shrnas(mrnas, params)
    results, _ = process_shrnas(cmap, mrnas, params)
    assert any(c.targets == 2 for c in results)
    assert any(c.targets == 1 for c in results)
    return results


async def test_target_all_filter(multi_mrna_candidates):
    """"Target all inputs" is on by default, restricting the view to
    candidates hitting every mRNA; 't' lifts that restriction, and
    toggling it again restores it."""
    app = ExploreApp(multi_mrna_candidates, total_mrnas=2, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        assert app._target_all_only is True
        n_restricted = len(app._filtered())
        assert n_restricted
        assert all(c.targets == 2 for c in app._filtered())
        await pilot.press("t")
        await pilot.pause()
        assert app._target_all_only is False
        filtered = app._filtered()
        assert len(filtered) > n_restricted
        # Toggling back on restores the restricted set.
        await pilot.press("t")
        await pilot.pause()
        assert app._target_all_only is True
        assert len(app._filtered()) == n_restricted


async def test_select_all_and_clear(candidates):
    """'a' selects every filtered candidate; shift+A clears the selection."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("a")
        await pilot.pause()
        assert len(app._selected) == len(app._filtered())
        await pilot.press("A")  # shift+a is uppercase
        await pilot.pause()
        assert len(app._selected) == 0


async def test_help_overlay_opens_and_closes(candidates):
    """'?' pushes the help screen; escape pops it back off the stack."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("question_mark")
        await pilot.pause()
        # Push happened: top screen is no longer the default ExploreApp screen.
        assert len(app.screen_stack) >= 2
        await pilot.press("escape")
        await pilot.pause()
        assert len(app.screen_stack) == 1


async def test_input_focus_suppresses_single_letter_binding(candidates):
    """Pressing 'm' while typing in the GC-min input must NOT toggle the motif
    filter (regression for the keybinding-collision class of bugs)."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        motif_before = app._motif_filter
        gc_min_input = app.query_one("#gc-min")
        gc_min_input.focus()
        await pilot.pause()
        await pilot.press("m")
        await pilot.pause()
        assert app._motif_filter == motif_before


async def test_oligo_screen_renders_for_each_backbone(candidates):
    """Every bundled backbone can be selected on the oligo screen and produces
    a non-empty preview, so none of them trip the config/preview pipeline."""
    from knockknock.backbone import load_backbones

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:3]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]
        for bb in load_backbones():
            oligo.query_one("#backbone-select").value = bb.name
            await pilot.pause()
            preview = oligo.query_one("#preview-text")
            rendered = str(preview.render())
            assert rendered.strip() != "", f"preview empty for {bb.name}"
            assert "Candidate 1/" in rendered


async def test_copy_action_falls_back_to_system_clipboard(candidates, monkeypatch):
    """`c` must reach a real clipboard even when OSC 52 is a no-op (e.g. macOS
    Terminal.app), so action_copy also shells out via _copy_to_system_clipboard."""
    from knockknock import oligo_tui as oligo_tui_mod

    calls: list[str] = []
    monkeypatch.setattr(oligo_tui_mod, "_copy_to_system_clipboard", calls.append)

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        await pilot.press("c")
        await pilot.pause()

        assert len(calls) == 1
        assert candidates[0].antisense in calls[0]


async def test_keep_sites_checkbox_reflects_and_edits_the_backbone_setting(
    candidates,
):
    """The checkbox is a stored backbone column, so it must load from the
    selected backbone and detach to "Custom" when the user changes it.

    pINDUCER ships with it on and pLKO with it off, so switching between them
    proves the box is driven by the database rather than left at whatever it
    happened to be — and that the programmatic sync does not itself trip
    `_mark_custom` (which is what the `prevent(..., Checkbox.Changed)` blocks
    are for).
    """
    from textual.widgets import Checkbox, Select

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test(size=(160, 50)) as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1], backbone_name="pINDUCER"))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        box = oligo.query_one("#keep-sites", Checkbox)
        assert box.value is True  # pINDUCER stores keep_sites=true
        assert oligo._selected_backbone_name == "pINDUCER"

        # Switching backbone re-seeds the box without marking the config custom.
        oligo.query_one("#backbone-select", Select).value = "pLKO"
        await pilot.pause()
        assert box.value is False
        assert oligo._selected_backbone_name == "pLKO"

        # A real edit does detach, and does reach the assembled oligos.
        box.value = True
        await pilot.pause()
        assert oligo._selected_backbone_name == "Custom"
        config = oligo._read_config()
        assert config.keep_sites is True
        assert config.site_bases(config.five_prime_overhang) == "T"  # AgeI
        assert config.site_bases(config.three_prime_overhang) == "G"  # EcoRI


async def test_copy_action_names_pair_like_fasta_export(candidates, monkeypatch):
    """The clipboard text must use the same name and layout as the FASTA
    export (`export_oligos_fasta_text`), not a reimplemented guide-prefix
    name — action_copy builds its text by calling that same formatter."""
    from knockknock import oligo_tui as oligo_tui_mod
    from knockknock.oligo import assemble_batch, export_oligos_fasta_text

    monkeypatch.setattr(oligo_tui_mod, "_copy_to_system_clipboard", lambda text: None)

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        copied: list[str] = []
        monkeypatch.setattr(oligo.app, "copy_to_clipboard", copied.append)

        oligo.action_copy()

        expected = export_oligos_fasta_text(
            assemble_batch(candidates[:1], oligo._config), oligo._config.backbone_name
        )
        assert copied == [expected]
        assert candidates[0].name in expected


async def test_save_backbone_flow(candidates):
    """Pressing `s`, naming the backbone, and confirming appends it to the CSV
    and makes it selectable without leaving the oligo screen."""
    from knockknock import backbone as backbone_mod

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        oligo.query_one("#flank-5").value = "AAAA"
        await pilot.pause()

        await pilot.press("s")
        await pilot.pause()
        modal = app.screen_stack[-1]
        assert modal is not oligo

        modal.query_one("#backbone-name-input").value = "Test Saved Backbone"
        await pilot.press("enter")
        await pilot.pause()

        assert app.screen_stack[-1] is oligo
        names = [b.name for b in backbone_mod.load_backbones()]
        assert "Test Saved Backbone" in names
        assert oligo.query_one("#backbone-select").value == "Test Saved Backbone"


async def test_backbone_button_switches_between_save_and_delete(candidates):
    """One always-visible button with two modes: Delete while a saved backbone
    is selected, Save once the fields diverge from it into "Custom". It is never
    hidden, so there is always a mouse-reachable way to manage the database."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]
        button = oligo.query_one("#save-backbone-btn")

        assert oligo.query_one("#backbone-select").value != "Custom"
        assert str(button.label) == "Delete Backbone"

        oligo.query_one("#flank-5").value = "AAAA"
        await pilot.pause()
        assert oligo.query_one("#backbone-select").value == "Custom"
        assert str(button.label) == "Save Backbone"


async def test_save_backbone_via_button_click(candidates):
    """The same save flow, but driven with a mouse click on the button
    instead of the `s` keybinding, for mouse-only accessibility."""
    from knockknock import backbone as backbone_mod

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        oligo.query_one("#flank-5").value = "AAAA"
        await pilot.pause()

        await pilot.click("#save-backbone-btn")
        await pilot.pause()
        modal = app.screen_stack[-1]
        assert modal is not oligo

        modal.query_one("#backbone-name-input").value = "Mouse Saved Backbone"
        await pilot.press("enter")
        await pilot.pause()

        assert app.screen_stack[-1] is oligo
        names = [b.name for b in backbone_mod.load_backbones()]
        assert "Mouse Saved Backbone" in names
        # Now that it is a saved backbone, the button offers to delete it.
        assert str(oligo.query_one("#save-backbone-btn").label) == "Delete Backbone"


async def test_save_backbone_rejects_invalid_config(candidates):
    """An invalid site/sequence must block the save, not silently write junk."""
    from knockknock import backbone as backbone_mod

    before_count = len(backbone_mod.load_backbones())

    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        oligo.query_one("#loop").value = "ACGN"  # degenerate code -> invalid
        await pilot.pause()

        await pilot.press("s")
        await pilot.pause()

        # No modal opened; still on the oligo screen with no new row saved.
        assert app.screen_stack[-1] is oligo
        assert len(backbone_mod.load_backbones()) == before_count


@pytest.mark.parametrize("widget_id", ["flank-5", "loop", "flank-3"])
@pytest.mark.parametrize(
    "bad_value,expected_snippet",
    [
        ("ACG1", "non-DNA character"),
        ("ACGN", "degenerate IUPAC code"),
    ],
)
async def test_non_dna_input_is_reported_not_crashed(
    candidates, widget_id, bad_value, expected_snippet
):
    """Typing a non-DNA or degenerate character must not crash the screen."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:1]))
        await pilot.pause()
        oligo: OligoScreen = app.screen_stack[-1]  # type: ignore[assignment]

        field = oligo.query_one(f"#{widget_id}")
        field.value = bad_value
        await pilot.pause()

        assert app.screen_stack[-1] is oligo  # still mounted, no crash
        preview = oligo.query_one("#preview-text")
        rendered = str(preview.render())
        assert expected_snippet in rendered
        assert "Candidate 1/" in rendered  # preview still renders the rest


# ---------------------------------------------------------------------------
#  Export format chooser
# ---------------------------------------------------------------------------


async def test_export_opens_format_chooser(candidates):
    """'e' pushes the format modal, offering one option per candidate format."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        assert isinstance(app.screen_stack[-1], ExportFormatModal)
        assert app.screen_stack[-1].query_one("#export-format-list").option_count == len(
            CANDIDATE_FORMATS
        )


async def test_export_chooser_escape_writes_nothing(candidates, tmp_path):
    """Cancelling the format modal returns to the main screen and leaves the
    output directory untouched — no partial or default-format file."""
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(tmp_path / "in.fa"))
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press("escape")
        await pilot.pause()
        assert app.screen_stack[-1] is app.screen
        assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "key,suffix",
    [("1", ".csv"), ("2", ".xlsx"), ("3", ".json"), ("4", ".fasta")],
)
async def test_export_writes_chosen_format(candidates, tmp_path, key, suffix):
    """Each number key in the chooser must produce a file of that format —
    beside the input FASTA, as resolve_export_path() promises."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(fasta))
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press(key)
        await pilot.pause()

    written = [p for p in tmp_path.iterdir() if p.suffix == suffix]
    assert len(written) == 1
    assert written[0].stat().st_size > 0


async def test_export_csv_keeps_targets_column(candidates, tmp_path):
    """The exported CSV includes a Targets column formatted as "n/total",
    not just the raw candidate count."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=2, input_path=str(fasta))
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press("1")
        await pilot.pause()

    csv_file = next(p for p in tmp_path.iterdir() if p.suffix == ".csv")
    content = csv_file.read_text()
    assert "Targets" in content
    assert "/2" in content


async def test_export_metadata_reflects_current_filters(candidates, tmp_path):
    """The provenance block records the filters in force in the TUI, not the
    wide-open ones the explore session started with."""
    import json as _json

    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(fasta))
    async with app.run_test() as pilot:
        await pilot.pause()
        app._gc_min = 25.0
        app._gc_max = 75.0
        app._motif_filter = False
        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press("3")  # JSON
        await pilot.pause()

    payload = _json.loads(
        next(p for p in tmp_path.iterdir() if p.suffix == ".json").read_text()
    )
    assert payload["meta"]["gc_min"] == 25.0
    assert payload["meta"]["gc_max"] == 75.0
    assert payload["meta"]["remove_motifs"] is False


async def test_rna_output_flag_reflected_in_table_and_export(candidates, tmp_path):
    """Regression: `--rna` used to be dropped entirely once `--tui`
    handed off to the TUI — neither the table nor its exports honored it."""
    import json as _json

    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(fasta), rna_output=True)
    async with app.run_test() as pilot:
        await pilot.pause()
        table = app.query_one("#results-table")
        row = table.get_row_at(0)
        guide_cell = str(row[2])
        assert "U" in guide_cell
        assert "T" not in guide_cell

        app.query_one("#results-table").focus()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press("3")  # JSON
        await pilot.pause()

    payload = _json.loads(
        next(p for p in tmp_path.iterdir() if p.suffix == ".json").read_text()
    )
    assert payload["meta"]["rna_output"] is True
    assert "U" in payload["candidates"][0]["guide_sequence"]


# ---------------------------------------------------------------------------
#  Top-N selection modal
# ---------------------------------------------------------------------------


async def test_select_top_n_flow(candidates):
    """'n', entering a count, and confirming selects that many top-ranked
    candidates from the current filtered/sorted view."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("n")
        await pilot.pause()

        assert len(app.screen_stack) >= 2
        modal = app.screen_stack[-1]
        modal.query_one("#top-n-input").value = "2"
        await pilot.press("enter")
        await pilot.pause()

        assert app.screen_stack[-1] is not modal
        assert len(app._selected) == 2
        top_2_names = {c.name for c in app._filtered()[:2]}
        assert app._selected == top_2_names


async def test_select_top_n_replaces_prior_selection(candidates):
    """Top-N *replaces* any pre-existing selection: previously checked rows
    that rank below the chosen N are unchecked, not kept as a union."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        # Pre-check two candidates that are NOT in the top 2 of the view.
        assert len(app._filtered()) >= 4
        app._selected = {c.name for c in app._filtered()[2:4]}
        app.query_one("#results-table").focus()
        await pilot.press("n")
        await pilot.pause()

        modal = app.screen_stack[-1]
        modal.query_one("#top-n-input").value = "2"
        await pilot.press("enter")
        await pilot.pause()

        top_2_names = {c.name for c in app._filtered()[:2]}
        assert app._selected == top_2_names


async def test_select_top_n_invalid_input_cancels(candidates):
    """A non-numeric top-N value dismisses the modal with no selection made,
    rather than crashing or selecting everything."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        app.query_one("#results-table").focus()
        await pilot.press("n")
        await pilot.pause()

        modal = app.screen_stack[-1]
        modal.query_one("#top-n-input").value = "not a number"
        await pilot.press("enter")
        await pilot.pause()

        assert app.screen_stack[-1] is app.screen
        assert len(app._selected) == 0


@pytest.mark.parametrize("key,suffix", [("1", ".fasta"), ("2", ".csv"), ("3", ".json")])
async def test_oligo_export_writes_candidates_and_oligos_files(
    candidates, tmp_path, key, suffix
):
    """Every format in the oligo screen's chooser writes two non-empty
    files of the matching type: one candidates file and one oligo pairs
    file — since a design done entirely from this screen shouldn't leave
    the guide-scoring data (score/dG/GC/targets) behind, and a single
    combined file would only work for some formats (CSV can't hold two
    tables)."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(fasta))
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(candidates[:2], input_path=str(fasta)))
        await pilot.pause()
        await pilot.press("e")
        await pilot.pause()
        assert isinstance(app.screen_stack[-1], ExportFormatModal)
        await pilot.press(key)
        await pilot.pause()

    written = [p for p in tmp_path.iterdir() if p.suffix == suffix]
    assert len(written) == 2
    assert all(p.stat().st_size > 0 for p in written)
    assert any("_candidates_" in p.name for p in written)
    assert any("_oligos_" in p.name for p in written)


@pytest.mark.parametrize("key,suffix", [("1", ".fasta"), ("2", ".csv"), ("3", ".json")])
async def test_oligo_export_candidates_file_matches_explorer_export(
    candidates, tmp_path, key, suffix
):
    """The oligo screen's candidates file must be produced the same way as
    the explorer screen's own `e` export (`export.write_candidates_export`),
    not a bespoke rendering — so it carries the same scoring fields
    (name/score/GC%/...) in the same shape."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(">x\nACGT\n")
    app = ExploreApp(candidates, total_mrnas=1, input_path=str(fasta))
    picked = candidates[:2]
    async with app.run_test() as pilot:
        await pilot.pause()
        app.push_screen(OligoScreen(picked, input_path=str(fasta)))
        await pilot.pause()
        await pilot.press("e")
        await pilot.pause()
        await pilot.press(key)
        await pilot.pause()

    candidates_out = next(
        p for p in tmp_path.iterdir() if p.suffix == suffix and "_candidates_" in p.name
    )
    text = candidates_out.read_text()
    for cand in picked:
        assert cand.name in text
    if suffix == ".fasta":
        for cand in picked:
            assert f">{cand.name} rank=" in text
            assert f"efficacy_score={cand.efficacy_score:.4f}" in text
    elif suffix == ".csv":
        assert "Guide_Sequence" in text
    elif suffix == ".json":
        assert '"candidates"' in text

    oligos_out = next(
        p for p in tmp_path.iterdir() if p.suffix == suffix and "_oligos_" in p.name
    )
    oligos_text = oligos_out.read_text()
    for cand in picked:
        assert cand.name in oligos_text


# ---------------------------------------------------------------------------
#  File picker: seeding, staying open on bad input, and the ctrl+o round trip
# ---------------------------------------------------------------------------


_GOOD_FASTA = ">x\n" + "ATGCGTAGCTAGCTAGCTAGCTGATCGTAGCTAGCTAGCT" + "\n"


async def _drive_picker(app: PickerApp, action) -> None:
    """Run ``action`` inside the picker and wait for the design worker only.

    Deliberately not ``app.workers.wait_for_complete()``: the `DirectoryTree`
    runs a loader worker that consumes a queue forever, so waiting on the
    whole group never returns. `WorkerCancelled` is expected rather than a
    failure — the happy path ends in ``app.exit()``, which cancels whatever is
    still registered.
    """
    async with app.run_test() as pilot:
        await pilot.pause()
        action(app)
        worker = next(w for w in app.workers if w.name == "_process_input")
        try:
            await worker.wait()
        except WorkerCancelled:
            pass
        await pilot.pause()


async def test_picker_designs_candidates_from_a_chosen_file(tmp_path):
    """Opening a valid FASTA runs the engine and stashes the results and the
    chosen path on the app, with `source_text` left unset (file, not paste)."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(_GOOD_FASTA)
    app = PickerApp()

    await _drive_picker(app, lambda a: a._open_file(fasta))

    assert app.candidates
    assert app.selected_path == fasta
    assert app.source_text is None


async def test_picker_stays_open_when_a_file_cannot_be_read(tmp_path):
    """A bad pick must not end the session — the next thing the user wants is
    to choose a different file, and after ctrl+o there is a loaded table to
    lose as well."""
    binary = tmp_path / "not_really.txt"
    binary.write_bytes(b"\xff\xfe\x00\x01\x02BINARY")
    app = PickerApp()

    await _drive_picker(app, lambda a: a._open_file(binary))

    assert app.candidates is None
    assert any("not a readable text file" in str(n.message) for n in app._notifications)


async def test_picker_stays_open_when_an_input_yields_no_candidates(tmp_path):
    """A readable FASTA that produces zero candidates (wildcards eating every
    window) leaves the picker up with a toast, not `app.candidates` set."""
    fasta = tmp_path / "wild.fa"
    fasta.write_text(">x\n" + "ATGCGTAGCN" * 5 + "\n")
    app = PickerApp()

    await _drive_picker(app, lambda a: a._open_file(fasta))

    assert app.candidates is None
    messages = [str(n.message) for n in app._notifications]
    assert any("No shRNA candidates" in m for m in messages)
    assert any("wildcard" in m for m in messages)


async def test_picker_seeds_the_tree_with_the_previously_opened_file(tmp_path):
    """The whole path stays on screen: the tree is rooted well above the file
    and expanded down to it, rather than rooted at the folder it sits in —
    which would hide everything around it and make the root row a bare,
    context-free folder name."""
    fasta = tmp_path / "nested" / "seed.fa"
    fasta.parent.mkdir()
    fasta.write_text(_GOOD_FASTA)
    app = PickerApp(initial_path=fasta.resolve())

    async with app.run_test() as pilot:
        for _ in range(3):
            await pilot.pause()
        tree = app.query_one("#picker-tree")

        root = Path(str(tree.path)).resolve()
        assert root == tree_root_for(fasta.resolve())
        assert root in fasta.resolve().parents
        assert tree.cursor_node is not None
        assert Path(tree.cursor_node.data.path).resolve() == fasta.resolve()

        # Every directory between the root and the file is expanded, so the
        # user can read the location off the tree and step back up it.
        shown = {
            Path(tree.get_node_at_line(line).data.path).resolve()
            for line in range(tree.last_line + 1)
        }
        assert fasta.parent.resolve() in shown
        assert tmp_path.resolve() in shown


async def test_directory_tree_go_up_reroots_to_the_parent_with_cursor_restored(
    tmp_path, monkeypatch
):
    """`Backspace` re-roots the tree one level up, with the folder just left
    under the cursor, expanded — the only way to reach anything above the
    picker's starting directory, since `DirectoryTree` can't navigate above
    its own root on its own, and collapsing that folder on the way up would
    hide the very listing the user was just looking at."""
    nested = tmp_path / "level1" / "level2"
    (nested / "child").mkdir(parents=True)
    monkeypatch.chdir(nested)
    app = PickerApp()

    async with app.run_test() as pilot:
        for _ in range(3):
            await pilot.pause()
        tree = app.query_one("#picker-tree")
        assert Path(tree.path).resolve() == nested.resolve()

        await pilot.press("backspace")
        for _ in range(3):
            await pilot.pause()
        assert Path(tree.path).resolve() == nested.parent.resolve()
        assert Path(tree.cursor_node.data.path).resolve() == nested.resolve()
        assert tree.cursor_node.is_expanded
        shown = {
            Path(tree.get_node_at_line(line).data.path).resolve()
            for line in range(tree.last_line + 1)
        }
        assert (nested / "child").resolve() in shown

        await pilot.press("backspace")
        for _ in range(3):
            await pilot.pause()
        assert Path(tree.path).resolve() == tmp_path.resolve()
        assert Path(tree.cursor_node.data.path).resolve() == nested.parent.resolve()


async def test_directory_tree_has_a_clickable_parent_row(tmp_path, monkeypatch):
    """The root's first child is a synthetic ".." row, right below the root's
    own label, so a mouse user can reach the parent folder without knowing
    about `Backspace`. Selecting it (which is what a click on the row does,
    same as `Enter`) re-roots the tree exactly like `action_go_up`."""
    nested = tmp_path / "level1" / "level2"
    nested.mkdir(parents=True)
    monkeypatch.chdir(nested)
    app = PickerApp()

    async with app.run_test() as pilot:
        for _ in range(3):
            await pilot.pause()
        tree = app.query_one("#picker-tree")
        parent_row = tree.root.children[0]
        assert parent_row.data.path == nested.parent
        assert not parent_row.allow_expand

        tree.cursor_line = parent_row.line
        await pilot.pause()
        await pilot.press("enter")
        for _ in range(3):
            await pilot.pause()
        assert Path(tree.path).resolve() == nested.parent.resolve()
        assert Path(tree.cursor_node.data.path).resolve() == nested.resolve()


async def test_directory_tree_cursor_starts_past_the_parent_row(tmp_path, monkeypatch):
    """On a fresh load with nothing to highlight, the cursor should land on
    real content, not on the ".." shortcut sitting above it."""
    nested = tmp_path / "level1" / "level2"
    (nested / "child").mkdir(parents=True)
    monkeypatch.chdir(nested)
    app = PickerApp()

    async with app.run_test() as pilot:
        for _ in range(3):
            await pilot.pause()
        tree = app.query_one("#picker-tree")
        assert not isinstance(tree.cursor_node.data, type(tree.root.children[0].data))
        assert tree.cursor_node.data.path == nested / "child"


async def test_enter_binding_to_open_is_visible_in_the_footer():
    """The picker's tree binds `Enter` to `select_cursor` with `show=True`,
    so it appears in the footer instead of staying hidden like `Tree`'s own
    default binding for the same key."""
    app = PickerApp()
    async with app.run_test() as pilot:
        for _ in range(3):
            await pilot.pause()
        tree = app.query_one("#picker-tree")
        tree.focus()
        await pilot.pause()
        enter_binding = app.screen.active_bindings["enter"].binding
        assert enter_binding.show is True
        assert enter_binding.action == "select_cursor"


async def test_tree_root_for_stays_at_home_for_a_file_inside_it():
    """A file under the home directory roots the tree at home, not its own
    parent folder — see `tree_root_for`'s docstring for why."""
    inside = Path.home() / "somewhere" / "targets.fa"
    assert tree_root_for(inside) == Path.home()


async def test_tree_root_for_falls_back_to_the_anchor_outside_home(tmp_path):
    """A file on an external volume or in /tmp is not under home; the nearest
    root that can still contain it is the filesystem anchor."""
    outside = tmp_path / "targets.fa"
    if Path.home() in outside.resolve().parents:
        pytest.skip("temp directory is inside home on this machine")
    assert tree_root_for(outside.resolve()) == Path(outside.resolve().anchor)


async def test_tree_root_for_defaults_to_cwd_without_a_file():
    """No file to seed with (e.g. first launch) roots the tree at the cwd,
    not home — `_FastaDirectoryTree.action_go_up` is what makes the rest of
    the filesystem, including home, reachable from there."""
    assert tree_root_for(None) == Path.cwd()


async def test_picker_seeds_the_paste_box_with_the_previously_pasted_text():
    """`initial_text` pre-fills the paste box, the paste half of re-seeding
    the picker after a round trip through `ctrl+o`."""
    app = PickerApp(initial_text=_GOOD_FASTA)

    async with app.run_test() as pilot:
        await pilot.pause()
        assert app.query_one("#paste-area").text == _GOOD_FASTA


async def test_picker_records_pasted_text_as_the_seed_for_next_time():
    """Submitting the paste box sets `source_text` (not `selected_path`) and
    `input_path` to "-", so `launch` knows to seed the next picker with text."""
    app = PickerApp()

    def paste_and_submit(a: PickerApp) -> None:
        a.query_one("#paste-area").text = _GOOD_FASTA
        a._submit_pasted_text()

    await _drive_picker(app, paste_and_submit)

    assert app.candidates
    assert app.source_text == _GOOD_FASTA
    assert app.selected_path is None
    assert app.input_path == "-"


async def test_ctrl_o_ends_the_explorer_with_an_open_request(candidates):
    """ctrl+o sets `open_requested` and exits the app, rather than pushing a
    screen — `ExploreApp` can't host `PickerApp`, so `launch`'s loop has to
    reopen the picker instead."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        assert app.open_requested is False
        await pilot.press("ctrl+o")
        await pilot.pause()
    assert app.open_requested is True


async def test_quitting_the_explorer_does_not_request_a_new_file(candidates):
    """'q' exits without setting `open_requested`, so `launch`'s loop ends
    the session instead of looping back to the picker."""
    app = ExploreApp(candidates, total_mrnas=1, input_path="test.fa")
    async with app.run_test() as pilot:
        await pilot.pause()
        await pilot.press("q")
        await pilot.pause()
    assert app.open_requested is False


async def test_explore_screen_shows_input_warnings_as_toasts(candidates):
    """`warnings` passed into `ExploreApp` (e.g. wildcard notes from loading)
    surface as notifications on mount, rather than being lost when the app
    takes over the screen."""
    app = ExploreApp(
        candidates,
        total_mrnas=1,
        input_path="test.fa",
        warnings=["Sequence 'x' contains 3 wildcard base(s) (N)."],
    )
    async with app.run_test() as pilot:
        await pilot.pause()
        assert any("wildcard base(s)" in str(n.message) for n in app._notifications)


async def test_launch_reopens_the_picker_seeded_with_the_current_file(
    monkeypatch, tmp_path
):
    """`launch` is what makes ctrl+o work: neither app can host the other, so
    the loop that runs them both is where "open another file" lives."""
    from knockknock import tui as tui_mod

    first = tmp_path / "first.fa"
    first.write_text(_GOOD_FASTA)
    second = tmp_path / "second.fa"
    second.write_text(_GOOD_FASTA)

    seeds: list[tuple] = []
    explored: list[str] = []

    class FakePicker:
        def __init__(self, initial_path=None, initial_text=None):
            seeds.append((initial_path, initial_text))
            self.candidates = ["candidate"]
            self.total_mrnas = 1
            self.input_path = str(second)
            self.warnings = []
            self.selected_path = second
            self.source_text = None

        def run(self):
            pass

    class FakeExplorer:
        runs = 0

        def __init__(self, candidates, total_mrnas, input_path, **kwargs):
            explored.append(input_path)
            FakeExplorer.runs += 1
            # Ask to reopen the first time round, then quit for real.
            self.open_requested = FakeExplorer.runs == 1

        def run(self):
            pass

    monkeypatch.setattr(tui_mod, "PickerApp", FakePicker)
    monkeypatch.setattr(tui_mod, "ExploreApp", FakeExplorer)

    tui_mod.launch(["candidate"], total_mrnas=1, input_path=str(first))

    assert explored == [str(first), str(second)]
    assert seeds == [(first.resolve(), None)]


async def test_launch_stops_when_the_picker_is_quit_without_a_choice(monkeypatch):
    """If the picker quits with `candidates` still `None`, `launch` returns
    immediately instead of opening the explorer with nothing to show."""
    from knockknock import tui as tui_mod

    class FakePicker:
        def __init__(self, initial_path=None, initial_text=None):
            self.candidates = None
            self.selected_path = None
            self.source_text = None

        def run(self):
            pass

    class FakeExplorer:
        def __init__(self, *args, **kwargs):
            raise AssertionError("explorer must not open without candidates")

    monkeypatch.setattr(tui_mod, "PickerApp", FakePicker)
    monkeypatch.setattr(tui_mod, "ExploreApp", FakeExplorer)

    tui_mod.launch()
