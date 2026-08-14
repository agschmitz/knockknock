"""Restriction cut sites and the single-stranded ends they leave behind.

Backbone junctions are specified the way cut sites are conventionally written:
the recognition sequence with a slash marking the **top-strand** cut, e.g.
``G/AATTC`` (EcoRI) or ``CTGCA/G`` (PstI).

That single slash is enough because Type II enzymes cut **symmetrically about
the centre of their recognition site**, so the bottom-strand cut is implied.
With site length ``L`` and ``c`` bases before the slash::

    k = L - 2c

    k > 0  ->  5' overhang of k nt,   top-strand sequence site[c : L-c]
    k = 0  ->  blunt
    k < 0  ->  3' overhang of |k| nt, top-strand sequence site[L-c : c]

Worked examples::

    G/AATTC     EcoRI    c=1  k=+4   5' AATT
    C/TCGAG     XhoI     c=1  k=+4   5' TCGA
    GG/CGCC     NarI     c=2  k=+2   5' CG
    /GATC       Sau3AI   c=0  k=+4   5' GATC
    CCC/GGG     SmaI     c=3  k= 0   blunt
    CTGCA/G     PstI     c=5  k=-4   3' TGCA
    GCATG/C     SphI     c=5  k=-4   3' CATG
    CATG/       NlaIII   c=4  k=-4   3' CATG

The ``CATG`` family shows why the slash position — not the recognition
sequence — is the load-bearing piece of information: ``/CATG`` leaves a 5'
``CATG``, ``C/ATG`` leaves a 5' ``AT``, and ``CATG/`` leaves a 3' ``CATG``.

**Everything here is stored in top-strand coordinates.** An `Overhang` records
the sequence of the single-stranded region as it reads 5'->3' on the top
strand, even when the protruding strand is physically the bottom one. That is
what removes the strand ambiguity that a bare "3' overhang = TCGA" string
carries, and it makes the drawing rules uniform: a top line drawn left-to-right
is the stored sequence verbatim, a bottom line drawn left-to-right is its
complement.

Enzymes that cut *outside* their recognition site (Type IIS: BsaI ``GGTCTC(1/5)``,
BbsI ``GAAGAC(2/6)``, BsmBI, SapI) cut asymmetrically and cannot be expressed
with a single slash — their overhang is whatever the plasmid happens to carry
in the spacer. For those, state the overhang directly: ``5'AATT``, ``3'TGCA``
or ``blunt``.

**A directly-stated overhang describes the digested backbone, not the insert**,
and — like any single strand — is written 5'->3' along the strand that actually
protrudes. Which strand that is depends on the junction::

    junction   kind   backbone strand that protrudes   stated 5'->3' runs
    left (5')  5'     bottom                           right-to-left
    left (5')  3'     top                              left-to-right
    right (3') 5'     top                              left-to-right
    right (3') 3'     bottom                           right-to-left

So converting a stated overhang into the top-strand coordinates stored on
`Overhang.top` is a reverse-complement in exactly the two rows where the
backbone's *bottom* strand protrudes, and the identity in the other two. The
insert then carries the complementary end, which is what makes ``5'TTAT`` and
``3'TTAT`` both put ``TTAT`` on the backbone rather than one of them landing on
the oligo. Real Type II overhangs are palindromic (their recognition sites are),
so this distinction only ever shows up on Type IIS / synthetic ends.

Cut-site notation is *not* transformed: ``G/AATTC`` names the top strand of the
uncut plasmid, which already is top-strand coordinates at either junction.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

FIVE_PRIME = "5p"
BLUNT = "blunt"
THREE_PRIME = "3p"

KINDS = (FIVE_PRIME, BLUNT, THREE_PRIME)

#: Which end of the insert a junction sits at. `LEFT` is the insert's 5'
#: junction (the ``5_prime_site`` column / ``--site-5`` flag), `RIGHT` its 3'
#: one. Needed wherever a *backbone-facing* overhang is read or written,
#: because the backbone protrudes with opposite strands at the two ends.
LEFT = "left"
RIGHT = "right"

JUNCTIONS = (LEFT, RIGHT)

#: Characters accepted as the cut marker in site notation.
SEPARATORS = "/^|"

#: The one DNA complement table in the codebase. Public because callers
#: outside this module need it for positional base-pair checks (the TUI's
#: duplex diagram) as well as for building complements — don't declare a
#: second copy elsewhere.
COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}

_EXPLICIT_RE = re.compile(r"^([53])['P]?-?([ACGT]+)$")

_BLUNT_WORDS = frozenset({"", "BLUNT", "NONE", "-"})

#: IUPAC ambiguity codes — a legitimate way to describe "one of several
#: bases" in a reference genome, but not something an oligo can be
#: synthesized from. Rejected wherever a synthesizable base is required
#: (overhangs, guides, flanks), but with a message distinct from other
#: non-DNA input (typos, whitespace, RNA's U, stray punctuation, ...).
#: Recognition sequences are the exception — see `validate_dna_sequence`.
DEGENERATE_BASES = frozenset("RYSWKMBDHVN")

#: Which concrete bases each IUPAC code stands for. Used to search a
#: sequence for a partly-degenerate recognition site (`site_pattern`); the
#: concrete bases map to themselves so one table covers both.
IUPAC_CLASSES: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def site_pattern(site: str) -> re.Pattern[str]:
    """Compile *site* into a regex that matches it base for base.

    Degenerate codes expand to a character class, so ``N/TCGAG``'s ``NTCGAG``
    matches ``CTCGAG`` and ``ATCGAG`` alike. A fully concrete site compiles to
    a plain literal, so this is a drop-in replacement for a substring search.
    """
    return re.compile("".join(f"[{IUPAC_CLASSES[b]}]" for b in site.upper()))


def site_occurs_in(site: str, sequence: str) -> bool:
    """True when *site* (possibly degenerate) occurs on either strand of
    *sequence*. `reverse_complement` is not degenerate-aware, so the reverse
    strand is covered by searching the sequence's own complement instead."""
    if not site:
        return False
    pattern = site_pattern(site)
    return bool(pattern.search(sequence)) or bool(
        pattern.search(reverse_complement(sequence))
    )


class CutSiteError(ValueError):
    """Raised when a cut-site / overhang specification cannot be parsed."""


def invalid_dna_chars(seq: str, *, allow_u: bool = False) -> tuple[str, str]:
    """Split the non-A/C/G/T characters of *seq* into (degenerate, other).

    Both are returned as sorted strings of the offending characters (each
    appearing once), so a caller can tell a degenerate IUPAC code apart from
    input that is not a DNA base at all. This is the single alphabet check
    shared by cut-site/overhang text, guide validation (`cli.py`), and mRNA
    FASTA validation (`cli.py`) — don't hand-roll another `set(...) -
    set("ACGT")` diff elsewhere; callers may still word their own error
    message differently for their audience.

    ``allow_u`` treats RNA's ``U`` as valid too, for callers validating raw
    (not-yet-DNA-normalized) input.
    """
    alphabet = set(COMPLEMENT) | ({"U"} if allow_u else set())
    bad = set(seq.upper()) - alphabet
    degenerate = "".join(sorted(bad & DEGENERATE_BASES))
    other = "".join(sorted(bad - DEGENERATE_BASES))
    return degenerate, other


def validate_dna_sequence(
    seq: str, label: str, *, allow_u: bool = False, allow_degenerate: bool = False
) -> None:
    """Raise `CutSiteError` unless *seq* is composed only of A/C/G/T (or also
    U, when ``allow_u`` is set).

    *other* (non-DNA, non-IUPAC) characters are reported ahead of degenerate
    codes when both are present, since they are the more surprising mistake.

    ``allow_degenerate`` permits IUPAC ambiguity codes. Recognition sequences
    set it: an enzyme whose site is only partly specified in the plasmid map
    (``N/TCGAG``) still pins the overhang down exactly, and the overhang is
    the part that has to be concrete. `Overhang.__post_init__` enforces that
    separately, on the overhang alone.
    """
    degenerate, other = invalid_dna_chars(seq, allow_u=allow_u)
    if other:
        raise CutSiteError(
            f"{label} contains non-DNA character(s) {other!r}. "
            "Only A, C, G, T are accepted."
        )
    if degenerate and not allow_degenerate:
        raise CutSiteError(
            f"{label} contains degenerate IUPAC code(s) {degenerate!r}. "
            "Give the concrete sequence — degenerate codes (N, R, Y, W, S, "
            "...) are not accepted here."
        )


def complement(seq: str) -> str:
    """Base-for-base complement, *without* reversing (positional complement)."""
    try:
        return "".join(COMPLEMENT[b] for b in seq)
    except KeyError as exc:  # pragma: no cover - guarded by parsing
        raise CutSiteError(f"Not a DNA base: {exc.args[0]!r}") from None


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


@dataclass(frozen=True)
class Overhang:
    """A single-stranded end, described in top-strand coordinates.

    ``kind`` is one of ``5p`` / ``blunt`` / ``3p`` and says *which strand
    protrudes* at this junction, given which end of the insert we are at:

    ========  ================================  ================================
    kind      at the insert's LEFT end           at the insert's RIGHT end
    ========  ================================  ================================
    ``5p``    top strand protrudes               bottom strand protrudes
    ``3p``    bottom strand protrudes            top strand protrudes
    ========  ================================  ================================

    ``top`` is always the sequence of that region read 5'->3' along the top
    strand, regardless of which strand is physically there.

    ``junction`` says which end of the insert this is (`LEFT` / `RIGHT`). It
    carries no weight for the insert-side geometry — `oligo.build_strands` and
    `oligo.make_duplex` already know which end they are placing — but it is
    what lets `backbone_top` state the same end the other way round: as the
    **backbone's** single-stranded overhang, written 5'->3'. That is the form
    the user types into a site field and the form `describe` shows back.
    """

    kind: str
    top: str = ""
    enzyme: str = ""
    #: Concrete recognition sequence, when the overhang came from site
    #: notation. Empty when the overhang was stated directly.
    recognition: str = ""
    #: The text this was parsed from, kept for round-tripping and display.
    spec: str = ""
    #: Which end of the insert this junction is — see `backbone_top`.
    junction: str = LEFT
    #: Recognition bases outside the overhang, on the left fragment's side of
    #: the cut (``C`` of ``C/TCGAG``). Empty for a directly stated overhang.
    site_head: str = ""
    #: Recognition bases outside the overhang, on the right fragment's side
    #: (``G`` of ``C/TCGAG``). Empty for a directly stated overhang.
    site_tail: str = ""

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise CutSiteError(
                f"Unknown overhang kind {self.kind!r}; expected one of {KINDS}."
            )
        if self.junction not in JUNCTIONS:
            raise CutSiteError(
                f"Unknown junction {self.junction!r}; expected one of {JUNCTIONS}."
            )
        if self.kind == BLUNT and self.top:
            raise CutSiteError("A blunt end cannot carry a sequence.")
        if self.kind != BLUNT and not self.top:
            raise CutSiteError(f"A {self.kind} overhang needs a sequence.")
        validate_dna_sequence(self.top, "Overhang sequence", allow_degenerate=True)
        degenerate, _ = invalid_dna_chars(self.top)
        if degenerate:
            # A degenerate code elsewhere in the recognition sequence is fine
            # (N/TCGAG pins the overhang down regardless); one *inside* the
            # overhang is not, because those are the bases that have to anneal.
            raise CutSiteError(
                f"The overhang would be {self.top!r}, which contains degenerate "
                f"IUPAC code(s) {degenerate!r}. The overhang itself must be "
                "concrete — the bases have to base-pair with the backbone. "
                "Degenerate codes are accepted elsewhere in a recognition "
                "sequence (e.g. 'N/TCGAG' still leaves a definite 5' TCGA)."
            )

    def __len__(self) -> int:
        return len(self.top)

    @property
    def is_blunt(self) -> bool:
        return self.kind == BLUNT

    @property
    def bottom(self) -> str:
        """The bottom strand of this region **as drawn left-to-right** (3'->5')."""
        return complement(self.top)

    @property
    def backbone_protrudes_on_top(self) -> bool:
        """True when it is the *backbone*'s top strand that sticks out here.

        The insert and the backbone protrude with opposite strands at the same
        junction, and which strand that is also flips between the two ends —
        so the backbone's top strand is the protruding one exactly when a 3'
        overhang sits at the left junction or a 5' overhang at the right one.
        """
        return (self.kind == THREE_PRIME) == (self.junction == LEFT)

    @property
    def backbone_top(self) -> str:
        """The backbone's single-stranded overhang, written 5'->3'.

        The user-facing form: what you read off a plasmid map, and what a site
        field is typed in. Equals `top` when the backbone's top strand is the
        protruding one, and its reverse complement when the bottom strand is —
        because a bottom strand runs 5'->3' right-to-left. `parse_overhang`
        applies the same transform in reverse to turn typed text into `top`.
        """
        if self.is_blunt:
            return ""
        return self.top if self.backbone_protrudes_on_top else reverse_complement(self.top)

    @property
    def insert_site_bases(self) -> str:
        """Recognition bases the *insert* must carry to rebuild the site.

        Digestion splits the recognition sequence across the two fragments.
        The overhang is shared, but the bases outside it are not: the left
        fragment keeps `site_head`, the right keeps `site_tail`. At the
        insert's left junction the insert *is* the right fragment, so it must
        supply `site_tail`; at its right junction it is the left fragment and
        must supply `site_head`.

        Ligation restores the site only if the insert brings these back —
        otherwise the junction is scarred and the enzyme will not re-cut. They
        are omitted by default (`OligoConfig.keep_sites`), since carrying them
        lengthens the insert and most users do not need to re-linearize.

        Empty when the junction was given as a bare overhang (``5'AATT``):
        there is no recognition sequence to rebuild. May contain a degenerate
        code, which means the insert's own contribution is undetermined — see
        `insert_site_bases_are_defined`.
        """
        if not self.recognition:
            return ""
        return self.site_tail if self.junction == LEFT else self.site_head

    @property
    def insert_site_bases_are_defined(self) -> bool:
        """True when `insert_site_bases` is concrete enough to synthesize.

        ``N/TCGAG`` is defined at a left junction (the insert supplies ``G``)
        but not at a right one (it would have to supply ``N``).
        """
        bases = self.insert_site_bases
        return bool(bases) and not any(b in DEGENERATE_BASES for b in bases)

    def describe(self) -> str:
        """Human-readable form, e.g. ``5' AATT`` / ``3' TGCA`` / ``blunt``.

        Shows the backbone-facing sequence (`backbone_top`), so it echoes what
        was typed into the site field rather than the internal top-strand
        coordinates — the two differ only for non-palindromic overhangs.
        """
        if self.is_blunt:
            return "blunt"
        return f"{self.kind[0]}' {self.backbone_top}"

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.describe()


@dataclass(frozen=True)
class CutSite:
    """A concrete recognition sequence plus its top-strand cut position."""

    site: str
    #: Number of bases before the slash on the top strand.
    top_cut: int
    enzyme: str = ""

    def __post_init__(self) -> None:
        if not self.site:
            raise CutSiteError("A cut site needs a recognition sequence.")
        validate_dna_sequence(
            self.site, "Recognition sequence", allow_degenerate=True
        )
        if not 0 <= self.top_cut <= len(self.site):
            raise CutSiteError(
                f"Cut position {self.top_cut} is outside the "
                f"{len(self.site)} bp recognition sequence."
            )

    @property
    def bottom_cut(self) -> int:
        """Bottom-strand cut, in top-strand coordinates (symmetric-cut rule)."""
        return len(self.site) - self.top_cut

    @property
    def protrusion(self) -> int:
        """``k = L - 2c``: >0 for 5' overhangs, 0 blunt, <0 for 3' overhangs."""
        return len(self.site) - 2 * self.top_cut

    @property
    def _bounds(self) -> tuple[int, int]:
        """The overhang's span in the site, as ``(start, end)``.

        The two cuts sit at `top_cut` and `bottom_cut`; whichever comes first
        opens the single-stranded region and the other closes it, so taking
        min/max covers 5' and 3' protrusions with one expression (and gives an
        empty span for a blunt cut, where the two coincide).
        """
        return min(self.top_cut, self.bottom_cut), max(self.top_cut, self.bottom_cut)

    @property
    def site_head(self) -> str:
        """Recognition bases before the overhang — kept by the left fragment."""
        return self.site[: self._bounds[0]]

    @property
    def site_tail(self) -> str:
        """Recognition bases after the overhang — kept by the right fragment."""
        return self.site[self._bounds[1] :]

    def overhang(self, junction: str = LEFT) -> Overhang:
        """The end this cut leaves, at *junction* of the insert.

        A method rather than a property because the junction has to come from
        the caller: the geometry below is junction-independent (the
        recognition sequence is given on the top strand either way), but the
        resulting `Overhang` needs to know its end to report a backbone-facing
        sequence.
        """
        k = self.protrusion
        start, end = self._bounds
        top = self.site[start:end]
        kind = FIVE_PRIME if k > 0 else THREE_PRIME if k < 0 else BLUNT
        return Overhang(
            kind=kind,
            top=top,
            enzyme=self.enzyme,
            recognition=self.site,
            spec=self.format(),
            junction=junction,
            site_head=self.site_head,
            site_tail=self.site_tail,
        )

    def format(self) -> str:
        """Render back to conventional notation, e.g. ``G/AATTC``."""
        return f"{self.site[: self.top_cut]}/{self.site[self.top_cut :]}"

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.format()


def _normalise(text: str) -> str:
    return "".join(text.split()).upper()


def parse_cut_site(text: str, enzyme: str = "") -> CutSite:
    """Parse conventional cut-site notation such as ``G/AATTC``.

    ``^`` and ``|`` are accepted as alternative cut markers. The marker may sit
    at either extreme (``/GATC``, ``CATG/``), but it must appear exactly once.
    """
    raw = _normalise(text)
    marks = [i for i, ch in enumerate(raw) if ch in SEPARATORS]
    if not marks:
        raise CutSiteError(
            f"No cut position in {text!r}. Mark the top-strand cut with a "
            "slash, e.g. 'G/AATTC'."
        )
    if len(marks) > 1:
        raise CutSiteError(
            f"Found {len(marks)} cut markers in {text!r}; expected exactly one. "
            "Only the top-strand cut is given — the bottom-strand cut is "
            "implied by the symmetry of the site."
        )
    cut = marks[0]
    site = raw[:cut] + raw[cut + 1 :]
    return CutSite(site=site, top_cut=cut, enzyme=enzyme)


def parse_overhang(text: str, enzyme: str = "", *, junction: str = LEFT) -> Overhang:
    """Parse a junction specification into an `Overhang`.

    Accepted forms::

        G/AATTC     cut-site notation; the overhang is derived
        5'AATT      a 5' overhang stated directly (Type IIS, synthetic ends)
        3'TGCA      a 3' overhang stated directly
        blunt       (also: an empty string)

    ``junction`` (`LEFT` / `RIGHT`) says which end of the insert this is, and
    is what makes a directly-stated overhang mean **the backbone's** end. The
    stated sequence is read 5'->3' along whichever backbone strand protrudes
    there, so it is reverse-complemented into `Overhang.top`'s top-strand
    coordinates whenever that strand is the bottom one — see the module
    docstring. Cut-site notation names the uncut plasmid's top strand and so
    is never transformed; ``junction`` is only recorded on the result.

    The default is `LEFT` so that a lone call reads naturally, but every
    caller handling a pair of junctions should pass both explicitly — passing
    `LEFT` for a 3' site silently designs an oligo that will not ligate.
    """
    raw = _normalise(text)
    if raw in _BLUNT_WORDS:
        return Overhang(kind=BLUNT, enzyme=enzyme, spec=text.strip(), junction=junction)

    if any(ch in raw for ch in SEPARATORS):
        return parse_cut_site(text, enzyme=enzyme).overhang(junction)

    match = _EXPLICIT_RE.match(raw)
    if match:
        end, seq = match.groups()
        kind = FIVE_PRIME if end == "5" else THREE_PRIME
        # The stated sequence runs 5'->3' along the backbone's protruding
        # strand; `top` is top-strand, left-to-right. Those agree only when
        # the backbone's top strand is the one sticking out.
        backbone_on_top = (kind == THREE_PRIME) == (junction == LEFT)
        return Overhang(
            kind=kind,
            top=seq if backbone_on_top else reverse_complement(seq),
            enzyme=enzyme,
            spec=text.strip(),
            junction=junction,
        )

    if raw and not set(raw) - set(COMPLEMENT):
        # A bare run of bases: almost certainly a recognition sequence whose
        # cut marker was left off. Without it the overhang is undetermined —
        # CATG alone could mean 5' CATG, 5' AT or 3' CATG.
        raise CutSiteError(
            f"No cut position in {text!r}. Mark the top-strand cut with a "
            f"slash — e.g. '{raw[:1]}/{raw[1:]}' — or state the overhang "
            f'directly as "5\'{raw}" or "3\'{raw}".'
        )

    raise CutSiteError(
        f"Could not read {text!r} as a junction. Give either cut-site notation "
        "with the top-strand cut marked (e.g. 'G/AATTC', 'CTGCA/G', 'CATG/'), "
        "or the overhang directly ('5'AATT', '3'TGCA', 'blunt')."
    )


def ends_are_interchangeable(left: Overhang, right: Overhang) -> bool:
    """True when the insert can also ligate in the reverse orientation.

    Flipping the insert 180 degrees makes its old bottom strand the new top
    strand, so the flipped left end presents ``reverse_complement(right.top)``.
    It still fits the backbone's left arm exactly when that equals
    ``left.top`` and the protrusion kinds match. Two blunt ends always match,
    which is why blunt cloning is never directional.

    Read backbone-side (`Overhang.backbone_top`) this is just "both junctions
    present the same overhang" — a vector cut with one enzyme at both ends
    reads identically at both, a directional EcoRI/XhoI pair does not.
    """
    if left.kind != right.kind:
        return False
    if left.is_blunt:
        return True
    return reverse_complement(right.top) == left.top
