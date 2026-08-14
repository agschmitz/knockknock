"""Tests for the shared sequence-input path used by every entry point."""

import pytest

from knockknock.seqinput import (
    SequenceInputError,
    load_mrnas,
    load_sequences,
    parse_sequence_text,
    read_sequence_file,
)


def test_parse_sequence_text_detects_fasta():
    """A '>'-led input is parsed as FASTA, id and sequence intact."""
    records = parse_sequence_text(">x\nATGCATGCATGCATGCATGCATGC\n")
    assert len(records) == 1
    assert records[0].id == "x"
    assert str(records[0].seq) == "ATGCATGCATGCATGCATGCATGC"


def test_parse_sequence_text_detects_multi_record_fasta():
    """Multiple '>'-headed records in one paste all come back, in order."""
    records = parse_sequence_text(">a\nATGC\n>b\nGGCC\n")
    assert [r.id for r in records] == ["a", "b"]


def test_parse_sequence_text_detects_strips_internal_whitespace_in_multi_record_fasta():
    """Multiple '>'-headed records in one paste all come back, in order."""
    records = parse_sequence_text(">a\nAT\nGC\n>b\nGG\nCC\n")
    assert [r.id for r in records] == ["a", "b"]
    assert str(records[0].seq) == "ATGC"
    assert str(records[1].seq) == "GGCC"


def test_parse_sequence_text_treats_headerless_input_as_one_raw_sequence():
    """No '>' header means the whole text is one sequence, id defaulted to
    'pasted_sequence'."""
    records = parse_sequence_text("ATGCATGCATGCATGCATGCATGC")
    assert len(records) == 1
    assert records[0].id == "pasted_sequence"
    assert str(records[0].seq) == "ATGCATGCATGCATGCATGCATGC"


def test_parse_sequence_text_strips_internal_whitespace_and_newlines():
    """A raw paste split across lines/spaces/tabs is joined into one sequence."""
    records = parse_sequence_text("ATGC ATGC\nATGC\tATGC\n")
    assert str(records[0].seq) == "ATGCATGCATGCATGC"


def test_parse_sequence_text_ignores_leading_blank_lines_before_fasta_header():
    """Blank lines before the '>' header don't defeat FASTA detection."""
    records = parse_sequence_text("\n\n>x\nATGCATGCATGCATGCATGCATGC\n")
    assert len(records) == 1
    assert records[0].id == "x"
    assert str(records[0].seq) == "ATGCATGCATGCATGCATGCATGC"


def test_parse_sequence_text_leading_blank_lines_without_header_still_raw():
    """Leading blank lines with no '>' anywhere still fall through to the
    raw-sequence path, not FASTA parsing."""
    records = parse_sequence_text("\n\nATGCATGCATGCATGCATGCATGC\n")
    assert len(records) == 1
    assert records[0].id == "pasted_sequence"
    assert str(records[0].seq) == "ATGCATGCATGCATGCATGCATGC"


def test_parse_sequence_text_treats_blank_input_as_no_records():
    """Blank text has no sequences at all — it is not one empty raw sequence,
    which would surface as a confusing "sequence 'pasted_sequence' is empty"."""
    assert parse_sequence_text("") == []
    assert parse_sequence_text("   \n\n\t ") == []


# ---------------------------------------------------------------------------
# load_mrnas / load_sequences — the one validation path every caller shares
# ---------------------------------------------------------------------------


def test_load_sequences_normalises_case_and_rna():
    """Lowercase RNA input comes back uppercased DNA — 'U' rewritten to 'T'
    — since the reference tables are DNA-based."""
    mrnas, warnings = load_sequences(">x\naugcaugcaugcaugcaugcaugc\n")
    assert mrnas == [("x", mrnas[0][1])]
    assert str(mrnas[0][1]) == "ATGCATGCATGCATGCATGCATGC"
    assert warnings == []


def test_load_sequences_rejects_duplicate_identifiers():
    """The engine keys candidates by identifier, so a duplicate would silently
    discard one input sequence entirely. This is the regression guard for the
    TUI path, which used to skip this check."""
    text = ">gene\nATGCATGCATGCATGCATGCATGC\n>gene\nGGCCGGCCGGCCGGCCGGCCGGCC\n"
    with pytest.raises(SequenceInputError, match="Duplicate FASTA identifier"):
        load_sequences(text)


def test_load_sequences_warns_about_wildcards_but_keeps_the_sequence():
    """Rejecting the record over an N run would throw away every well-defined
    region with it; the wildcards ride through and cost only the guide windows
    that overlap them."""
    mrnas, warnings = load_sequences(">x\nATGCNATGCNATGCRATGCN\n")

    assert str(mrnas[0][1]) == "ATGCNATGCNATGCRATGCN"
    assert len(warnings) == 1
    assert "4 wildcard base(s)" in warnings[0]
    assert "N, R" in warnings[0]


def test_load_sequences_warns_per_record_only_where_wildcards_occur():
    """With multiple records, only the one carrying a wildcard gets a warning
    — a clean sibling record doesn't."""
    text = ">clean\nATGCATGCATGCATGCATGCATGC\n>dirty\nATGCATGCNTGCATGCATGCATGC\n"
    _, warnings = load_sequences(text)

    assert len(warnings) == 1
    assert "'dirty'" in warnings[0]


def test_load_sequences_still_rejects_non_nucleotide_characters():
    """Wildcards are meaningful; a stray '@' is a broken file."""
    with pytest.raises(SequenceInputError, match="invalid characters"):
        load_sequences(">x\nATGC@ATGC@ATGC@ATGC@\n")


def test_load_sequences_rejects_empty_record():
    """A record with a header but no sequence is rejected, even when other
    records in the same input are fine."""
    with pytest.raises(SequenceInputError, match="is empty"):
        load_sequences(">x\n\n>y\nATGCATGCATGCATGCATGCATGC\n")


def test_load_sequences_rejects_blank_input():
    """Empty text yields no records, which load_mrnas reports as such rather
    than as an empty-sequence error."""
    with pytest.raises(SequenceInputError, match="No sequences found"):
        load_sequences("")


def test_load_mrnas_rejects_empty_record_list():
    """load_mrnas itself, not just load_sequences, rejects an empty record list."""
    with pytest.raises(SequenceInputError, match="No sequences found"):
        load_mrnas([])


def test_read_sequence_file_reports_binary_input_as_an_input_error(tmp_path):
    """A binary file is a realistic pick in the TUI's file browser, where
    ".txt" is an accepted suffix — it must not escape as a UnicodeDecodeError."""
    binary = tmp_path / "not_really.txt"
    binary.write_bytes(b"\xff\xfe\x00\x01\x02BINARY")

    with pytest.raises(SequenceInputError, match="not a readable text file"):
        read_sequence_file(binary)


def test_read_sequence_file_reports_missing_file_as_an_input_error(tmp_path):
    """A nonexistent path surfaces as SequenceInputError, not a raw OSError."""
    with pytest.raises(SequenceInputError, match="Could not read"):
        read_sequence_file(tmp_path / "nope.fa")
