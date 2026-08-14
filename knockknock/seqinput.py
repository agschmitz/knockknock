"""Sequence input — the one path from raw text to engine-ready mRNAs.

Every way a sequence can enter the tool (a FASTA file on the command line,
stdin, the TUI's file picker, the TUI's paste box, a file dragged onto the
terminal) funnels through `load_sequences` here. There is deliberately no
second implementation: the CLI and the TUI used to parse and validate
separately, which is how the TUI ended up silently dropping records whose
FASTA header was a duplicate of an earlier one.

The pipeline is:

    text -> parse_sequence_text()  -> [SeqRecord]        # FASTA or one raw sequence
         -> load_mrnas()         -> LoadedSequences    # validated, DNA-normalised

`load_sequences` is the two of them together, and is what callers should use.
Anything rejected raises `SequenceInputError` with a message written for an
end user; callers render it however suits their interface. Anything merely
*worth mentioning* comes back in `LoadedSequences.warnings` instead — see
that class for why wildcard bases are in the second category.
"""

import re
from io import StringIO
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from knockknock.cutsite import invalid_dna_chars


class SequenceInputError(ValueError):
    """Raised when input text cannot be turned into usable mRNA records."""


class LoadedSequences(NamedTuple):
    """Engine-ready mRNAs, plus non-fatal notes about how they were read.

    A pair rather than a bare list of mRNAs because degenerate IUPAC codes
    (``N``, ``R``, ``Y``, ...) are **not** fatal: the sequence is kept
    verbatim, and the guide windows that overlap a wildcard are dropped later
    by `engine.score_one`, so the defined A/C/G/T stretches around them can
    still be targeted. That trade-off is invisible unless someone says so, so
    every caller is handed a warning to render. Unpack it::

        mrnas, warnings = load_sequences(text)
    """

    mrnas: list[tuple[str, Seq]]
    warnings: list[str]


def parse_sequence_text(text: str) -> list[SeqRecord]:
    """Parse ``text`` as FASTA if it has a '>' header line, else as one raw,
    unheadered sequence (whitespace and line breaks stripped).

    Blank or whitespace-only text yields no records at all, rather than one
    empty raw sequence — so an empty file reports "no sequences found" instead
    of an "empty sequence" named after the raw-paste placeholder.

    Raises whatever ``Bio.SeqIO`` raises if a '>'-led input isn't valid FASTA.
    """
    text = text.lstrip("\n\r \t")
    if not text.strip():
        return []

    if re.search(r"^\s*>", text, re.MULTILINE):
        return list(SeqIO.parse(StringIO(text), "fasta"))

    seq_str = re.sub(r"\s+", "", text)
    return [SeqRecord(Seq(seq_str), id="pasted_sequence", description="")]


def load_mrnas(records: list[SeqRecord]) -> LoadedSequences:
    """Validate parsed records and normalise them into engine input.

    Returns `LoadedSequences`: ``(id, Seq)`` pairs with every sequence
    uppercased and RNA's ``U`` rewritten to ``T`` — the scoring reference
    tables are DNA-based, so this is the one place that normalisation happens
    — plus a warning per record that carries wildcard bases.

    Raises `SequenceInputError` when a record is empty, carries a character
    that is neither a nucleotide nor an IUPAC ambiguity code, or repeats an
    earlier record's identifier. Duplicate identifiers are fatal rather than a
    warning because the engine keys its candidate map by identifier: a
    silently-merged pair of records loses one sequence's candidates entirely
    while still reporting two inputs.

    **Degenerate IUPAC codes are kept, not rejected.** A reference or assembly
    sequence routinely carries an ``N`` run, and rejecting the whole record
    over it throws away every well-defined region as well. The wildcards ride
    through to `engine.make_shrnas`, whose windows are then dropped one by one
    in `engine.score_one` if they overlap one — so only the guides that
    actually touch an undefined base are lost.
    """
    if not records:
        raise SequenceInputError("No sequences found in input.")

    seen_ids: set[str] = set()
    mrnas: list[tuple[str, Seq]] = []
    warnings: list[str] = []
    for record in records:
        seq_str = str(record.seq).upper()
        if not seq_str:
            raise SequenceInputError(f"Sequence {record.id!r} is empty.")

        degenerate, other = invalid_dna_chars(seq_str, allow_u=True)
        if other:
            raise SequenceInputError(
                f"Sequence {record.id!r} contains invalid characters: "
                f"{other!r}. Expected nucleotides A, C, G, T, or U, or an "
                "IUPAC ambiguity code (N, R, Y, W, S, ...)."
            )

        if record.id in seen_ids:
            raise SequenceInputError(
                f"Duplicate FASTA identifier {record.id!r}. Give every "
                "sequence a unique header."
            )
        seen_ids.add(record.id)

        if degenerate:
            count = sum(seq_str.count(code) for code in degenerate)
            warnings.append(
                f"Sequence {record.id!r} contains {count} wildcard base(s) "
                f"({', '.join(degenerate)}). Guides overlapping them are "
                "discarded — only the defined A/C/G/T stretches are targeted."
            )

        mrnas.append((record.id, Seq(seq_str.replace("U", "T"))))

    return LoadedSequences(mrnas, warnings)


def load_sequences(text: str) -> LoadedSequences:
    """Parse and validate ``text`` into engine-ready mRNAs.

    The single entry point for sequence input. ``text`` is either a FASTA
    (one or more records) or a bare unheadered sequence — see
    `parse_sequence_text` for how the two are told apart.
    """
    try:
        records = parse_sequence_text(text)
    except SequenceInputError:
        raise
    except Exception as exc:
        raise SequenceInputError(f"Could not parse input as FASTA: {exc}") from exc
    return load_mrnas(records)


def read_sequence_file(path: Path) -> str:
    """Read a sequence file as text, with errors phrased for an end user.

    Binary or non-UTF-8 files are a realistic pick in the TUI's file browser
    (``.txt`` is an accepted suffix), so the decode failure is turned into a
    `SequenceInputError` rather than an unhandled `UnicodeDecodeError`.
    """
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise SequenceInputError(
            f"{path} is not a readable text file — it looks binary, or is not "
            "UTF-8 encoded. Give a plain-text FASTA file."
        ) from exc
    except OSError as exc:
        raise SequenceInputError(f"Could not read {path}: {exc}") from exc
