"""Tests for the backbone CSV loader and the save-back-to-CSV path."""

import pytest

from knockknock import backbone as backbone_mod
from knockknock.backbone import (
    Backbone,
    backbone_label,
    delete_backbone,
    get_backbone,
    load_backbones,
    restore_bundled_backbones,
    save_backbone,
)
from knockknock.cutsite import CutSiteError, parse_overhang
from knockknock.oligo import OligoConfig


@pytest.fixture
def isolated_backbones_csv(sandbox_backbone_database):
    """The throwaway backbone database for this test.

    The isolation itself is done by the autouse `sandbox_backbone_database`
    fixture in conftest.py, which every test gets. This is just a readable alias
    for tests that want to inspect or edit the CSV directly.
    """
    return sandbox_backbone_database


def _custom_backbone(name: str) -> Backbone:
    return Backbone(
        name=name,
        five_prime_enzyme="",
        five_prime_overhang=parse_overhang("G/AATTC"),
        five_prime_flank="AAA",
        loop="TAGTGAAGCCACAGATGTA",
        three_prime_flank="CCC",
        three_prime_enzyme="",
        three_prime_overhang=parse_overhang("C/TCGAG"),
        guide_on_5p=False,
    )


def test_save_backbone_appends_row_and_refreshes_cache(isolated_backbones_csv):
    """A saved backbone is immediately visible through `load_backbones` — the
    in-memory cache is invalidated, not left stale after the write."""
    before = len(load_backbones())

    save_backbone(_custom_backbone("My Custom Backbone"))

    after = load_backbones()
    assert len(after) == before + 1
    saved = next(b for b in after if b.name == "My Custom Backbone")
    assert saved.five_prime_flank == "AAA"
    assert saved.loop == "TAGTGAAGCCACAGATGTA"
    assert saved.three_prime_flank == "CCC"
    assert saved.guide_on_5p is False


def test_save_backbone_round_trips_through_disk(isolated_backbones_csv):
    """Bypassing the cache and re-reading the CSV from scratch still recovers
    the overhang specs exactly, proving the write itself (not just the cache
    update) is correct."""
    save_backbone(_custom_backbone("Round Trip Backbone"))

    # Force a fresh read from disk rather than the in-memory cache.
    backbone_mod._BACKBONES = None
    reloaded = load_backbones()
    saved = next(b for b in reloaded if b.name == "Round Trip Backbone")
    assert saved.five_prime_overhang.spec == "G/AATTC"
    assert saved.three_prime_overhang.spec == "C/TCGAG"


def test_keep_sites_round_trips_through_the_csv(isolated_backbones_csv):
    """`keep_sites` is a stored column, so it must survive a write/re-read."""
    backbone = _custom_backbone("Keeps Its Sites")
    backbone.keep_sites = True
    save_backbone(backbone)

    backbone_mod._BACKBONES = None
    saved = next(b for b in load_backbones() if b.name == "Keeps Its Sites")
    assert saved.keep_sites is True
    assert OligoConfig.from_backbone(saved).keep_sites is True


def test_bundled_keep_sites_defaults(isolated_backbones_csv):
    """pINDUCER regenerates its XhoI/EcoRI sites; the rest do not by default.

    pINDUCER's 5' flank deliberately no longer starts with the `G` XhoI needs
    — that base is now supplied by the site restoration instead of being
    duplicated in the flank, so the two cannot disagree.
    """
    assert get_backbone("pINDUCER").keep_sites is True
    assert not get_backbone("pINDUCER").five_prime_flank.startswith("G")
    for name in ("pGIPZ", "pTRIPZ", "pLKO", "pHREP"):
        assert get_backbone(name).keep_sites is False, name


def test_a_database_without_the_keep_sites_column_still_loads(
    isolated_backbones_csv, tmp_path, monkeypatch
):
    """A backbones.csv written before the column existed must not KeyError.

    The user's own database under ~/.config is seeded once and then only
    appended to, so an installed copy can predate any column added later.
    """
    legacy = tmp_path / "legacy.csv"
    legacy.write_text(
        "name,5_prime_enzyme,5_prime_site,5_prime_flank,loop,3_prime_flank,"
        "3_prime_enzyme,3_prime_site,guide_on_5p\n"
        "Legacy,XhoI,C/TCGAG,AAA,TTT,CCC,EcoRI,G/AATTC,false\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("KNOCKKNOCK_BACKBONES_PATH", str(legacy))
    backbone_mod._BACKBONES = None

    loaded = load_backbones()
    assert [b.name for b in loaded] == ["Legacy"]
    assert loaded[0].keep_sites is False
    assert loaded[0].mismatch_positions == [1]


def test_mismatch_positions_round_trips_through_the_csv(isolated_backbones_csv):
    """`mismatch_positions` is a stored column, so it must survive a
    write/re-read, the same way `keep_sites` does."""
    backbone = _custom_backbone("Custom Mismatch Backbone")
    backbone.mismatch_positions = [2, 4]
    save_backbone(backbone)

    backbone_mod._BACKBONES = None
    saved = next(
        b for b in load_backbones() if b.name == "Custom Mismatch Backbone"
    )
    assert saved.mismatch_positions == [2, 4]
    assert OligoConfig.from_backbone(saved).mismatch_positions == [2, 4]


def test_bundled_mismatch_position_defaults(isolated_backbones_csv):
    """Every bundled backbone stores the conventional position 1 — this
    feature changes nothing about existing behaviour until a backbone
    actually declares something else."""
    for name in ("pINDUCER", "pGIPZ", "pTRIPZ", "Open Biosystems", "pLKO", "pHREP"):
        assert get_backbone(name).mismatch_positions == [1], name


def test_a_database_without_the_mismatch_positions_column_still_loads(
    isolated_backbones_csv, tmp_path, monkeypatch
):
    """A backbones.csv written before this column existed must not KeyError,
    the same guarantee `keep_sites` already has."""
    legacy = tmp_path / "legacy.csv"
    legacy.write_text(
        "name,5_prime_enzyme,5_prime_site,5_prime_flank,loop,3_prime_flank,"
        "3_prime_enzyme,3_prime_site,guide_on_5p,keep_sites\n"
        "Legacy,XhoI,C/TCGAG,AAA,TTT,CCC,EcoRI,G/AATTC,false,false\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("KNOCKKNOCK_BACKBONES_PATH", str(legacy))
    backbone_mod._BACKBONES = None

    loaded = load_backbones()
    assert [b.name for b in loaded] == ["Legacy"]
    assert loaded[0].mismatch_positions == [1]


def test_mismatch_positions_zero_disables(isolated_backbones_csv, tmp_path, monkeypatch):
    """A backbone row can store `0` (no mismatch) — same sentinel convention
    as `-m 0`/`--mismatch 0`."""
    csv_path = tmp_path / "backbones.csv"
    csv_path.write_text(
        "name,5_prime_enzyme,5_prime_site,5_prime_flank,loop,3_prime_flank,"
        "3_prime_enzyme,3_prime_site,guide_on_5p,keep_sites,mismatch_positions\n"
        "NoMismatch,XhoI,C/TCGAG,AAA,TTT,CCC,EcoRI,G/AATTC,false,false,0\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("KNOCKKNOCK_BACKBONES_PATH", str(csv_path))
    backbone_mod._BACKBONES = None

    assert get_backbone("NoMismatch").mismatch_positions == []


def test_load_backbones_reports_line_number_on_bad_mismatch_positions(
    isolated_backbones_csv,
):
    """A malformed `mismatch_positions` cell raises with the CSV line number
    and backbone name, the same as a malformed junction spec does."""
    with open(isolated_backbones_csv, "a", newline="") as f:
        f.write(
            "Bad Mismatch,,G/AATTC,AAA,TAGTGAAGCCACAGATGTA,CCC,,C/TCGAG,"
            "false,false,abc\n"
        )

    with pytest.raises((CutSiteError, ValueError)) as excinfo:
        load_backbones()

    message = str(excinfo.value)
    assert "Bad Mismatch" in message
    assert "line" in message


def test_save_backbone_rejects_duplicate_name(isolated_backbones_csv):
    """Saving under a name that already exists raises instead of silently
    creating a second row with that name; callers must `delete_backbone` first
    to replace one."""
    existing_name = load_backbones()[0].name

    with pytest.raises(ValueError, match="already exists"):
        save_backbone(_custom_backbone(existing_name))


def test_save_backbone_seeds_a_configured_user_database(tmp_path, monkeypatch):
    """When `KNOCKKNOCK_BACKBONES_PATH` points at a file that doesn't exist yet,
    the first save creates it seeded with the bundled defaults, not just the
    one new row — so a fresh configured database still has the shipped
    backbones alongside the user's own."""
    user_csv = tmp_path / "config" / "backbones.csv"
    monkeypatch.setenv("KNOCKKNOCK_BACKBONES_PATH", str(user_csv))
    monkeypatch.setattr(backbone_mod, "_BACKBONES", None)

    save_backbone(_custom_backbone("Configured User Backbone"))

    assert user_csv.exists()
    assert "pINDUCER" in user_csv.read_text()
    assert "Configured User Backbone" in user_csv.read_text()


def test_backbone_label_appends_differing_enzymes(isolated_backbones_csv):
    """The bare CSV name gets its two junction enzymes appended, 5' then 3'."""
    bb = get_backbone("pINDUCER")
    assert backbone_label(bb) == "pINDUCER (XhoI/EcoRI)"


def test_backbone_label_collapses_a_shared_enzyme(isolated_backbones_csv):
    """A backbone cut by the same enzyme at both junctions (e.g. a Golden Gate
    design) shows it once, not as "Esp3I/Esp3I"."""
    bb = get_backbone("pHREP")
    assert backbone_label(bb) == "pHREP (Esp3I)"


def test_backbone_label_omits_parens_without_enzymes(isolated_backbones_csv):
    """A custom backbone saved with no enzyme recorded (e.g. from the TUI's
    'Save Backbone' flow) shows just the bare name — no empty "()"."""
    bb = _custom_backbone("My Custom Backbone")
    assert backbone_label(bb) == "My Custom Backbone"


def test_load_backbones_reports_line_number_on_bad_site(isolated_backbones_csv):
    """A malformed junction spec should raise with the CSV line number and
    backbone name, not a bare parser error, so a hand-edited CSV is easy to
    fix."""
    with open(isolated_backbones_csv, "a", newline="") as f:
        # No slash in the site: parse_overhang can't derive an overhang.
        f.write("Bad Backbone,,GAATTC,AAA,TAGTGAAGCCACAGATGTA,CCC,,C/TCGAG,false\n")

    with pytest.raises(CutSiteError) as excinfo:
        load_backbones()

    message = str(excinfo.value)
    assert "Bad Backbone" in message
    assert "line" in message


# ---------------------------------------------------------------------------
# delete_backbone / restore_bundled_backbones
# ---------------------------------------------------------------------------


def test_delete_backbone_removes_a_saved_backbone(isolated_backbones_csv):
    """A backbone the user saved can be removed again from the database."""
    save_backbone(_custom_backbone("Doomed Backbone"))
    assert "Doomed Backbone" in [b.name for b in load_backbones()]

    delete_backbone("Doomed Backbone")

    assert "Doomed Backbone" not in [b.name for b in load_backbones()]


def test_delete_backbone_also_removes_bundled_ones(isolated_backbones_csv):
    """Any backbone can be deleted, bundled defaults included — the user's own
    database is theirs to curate."""
    bundled = next(iter(backbone_mod._bundled_backbone_names()))

    delete_backbone(bundled)

    assert bundled not in [b.name for b in load_backbones()]


def test_delete_backbone_leaves_the_other_rows_intact(isolated_backbones_csv):
    """Deleting rewrites the whole CSV, so the untouched rows must survive
    byte-for-byte in their parsed form."""
    before = [b.name for b in load_backbones()]
    victim = before[0]

    delete_backbone(victim)

    assert [b.name for b in load_backbones()] == before[1:]


def test_delete_backbone_rejects_an_unknown_name(isolated_backbones_csv):
    """Deleting something that isn't there is an error, not a silent no-op."""
    with pytest.raises(ValueError, match="No backbone named"):
        delete_backbone("Never Existed")


def test_delete_backbone_never_touches_the_bundled_package_file(isolated_backbones_csv):
    """Deletion must edit the user's copy only. The file inside the installed
    package is shared, often read-only, and is what makes restore possible."""
    bundled_path = backbone_mod._bundled_data_path()
    before = bundled_path.read_text(encoding="utf-8")

    delete_backbone(next(iter(backbone_mod._bundled_backbone_names())))

    assert bundled_path.read_text(encoding="utf-8") == before


def test_restore_bundled_backbones_brings_deleted_defaults_back(isolated_backbones_csv):
    """The counterpart to unrestricted deletion: shipped definitions are never
    truly lost."""
    bundled = sorted(backbone_mod._bundled_backbone_names())[:2]
    for name in bundled:
        delete_backbone(name)
    assert not set(bundled) & {b.name for b in load_backbones()}

    restored = restore_bundled_backbones()

    assert sorted(restored) == bundled
    assert set(bundled) <= {b.name for b in load_backbones()}


def test_restore_bundled_backbones_keeps_the_users_own_backbones(isolated_backbones_csv):
    """Restoring must not wipe custom work — it only re-adds what is missing."""
    save_backbone(_custom_backbone("My Own Backbone"))
    delete_backbone(next(iter(backbone_mod._bundled_backbone_names())))

    restore_bundled_backbones()

    assert "My Own Backbone" in [b.name for b in load_backbones()]


def test_restore_bundled_backbones_is_a_no_op_when_nothing_is_missing(
    isolated_backbones_csv,
):
    """Running restore on an untouched database changes nothing."""
    before = [b.name for b in load_backbones()]

    assert restore_bundled_backbones() == []
    assert [b.name for b in load_backbones()] == before
