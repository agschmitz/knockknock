"""Shared fixtures.

The autouse fixtures here sandbox everything the tool can write to, so no test
can reach the developer's real home directory or leave files in the repository.
That matters more than usual for this package: exports default to the current
working directory, and saved backbones default to ``~/.config/knockknock/``.
"""

import shutil

import pytest

from knockknock import backbone as backbone_mod
from knockknock import tui as tui_mod


@pytest.fixture(autouse=True)
def sandbox_working_directory(tmp_path, monkeypatch):
    """Run every test in its own empty directory.

    Exports with no explicit ``-o`` land in the current working directory, so
    without this a test run would scatter timestamped CSV/XLSX files through the
    repository. It also means a test can assert on ``tmp_path.iterdir()`` to
    find what an export wrote.
    """
    monkeypatch.chdir(tmp_path)


@pytest.fixture(autouse=True)
def sandbox_backbone_database(tmp_path_factory, monkeypatch):
    """Point the backbone database at a throwaway copy of the bundled CSV.

    `save_backbone`/`delete_backbone` write to ``~/.config/knockknock/`` by
    default. Without this, running the suite would add test backbones to the
    developer's own database — and the tests would then see different data on
    the next run. The module-level cache is cleared either side so the copy is
    actually read.

    Deliberately NOT inside ``tmp_path``: that is the working directory tests
    inspect to see what an export wrote, and a stray database folder there would
    show up as an unexpected file.
    """
    db = tmp_path_factory.mktemp("backbone_db") / "backbones.csv"
    shutil.copyfile(backbone_mod._bundled_data_path(), db)
    monkeypatch.setenv("KNOCKKNOCK_BACKBONES_PATH", str(db))
    monkeypatch.setattr(backbone_mod, "_BACKBONES", None)
    yield db
    backbone_mod._BACKBONES = None


@pytest.fixture(autouse=True)
def sandbox_tui_state(tmp_path_factory, monkeypatch):
    """Keep the TUI's onboarding-banner state out of the real config directory.

    Outside ``tmp_path`` for the same reason as the backbone database above.
    """
    monkeypatch.setattr(
        tui_mod, "_STATE_FILE", tmp_path_factory.mktemp("tui_state") / "state.json"
    )
