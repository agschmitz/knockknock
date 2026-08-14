"""KnockKnock — Design shRNA guide sequences for miR30-based knockdown."""

from importlib.metadata import PackageNotFoundError, version

try:
    # Single source of truth: the version declared in pyproject.toml, read back
    # from the installed distribution's metadata. Previously this was a literal
    # here *and* in pyproject.toml, which is two places to forget on a release.
    __version__ = version("knockknock-shrna")
except PackageNotFoundError:  # pragma: no cover - running from an uninstalled tree
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
