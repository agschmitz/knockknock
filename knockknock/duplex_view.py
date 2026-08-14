"""ASCII-art rendering of the insert duplex sitting in its digested backbone.

Pure presentation code: every function here is a function of its arguments
only. Nothing in this module touches a widget or imports Textual, so the
diagram can be rendered (and tested) outside the TUI — and reused by any other
front end that can display colour-tagged text.

Colours are emitted as Rich console markup (``[#rrggbb]...[/#rrggbb]``), which
is what `oligo_tui.OligoScreen` feeds to a `Static`. The plain-text geometry
(column offsets, padding, label placement) is independent of that markup.
"""

from knockknock.cutsite import (
    COMPLEMENT,
    FIVE_PRIME,
    THREE_PRIME,
    Overhang,
    complement,
)
from knockknock.oligo import OligoConfig, make_duplex

# Muted terminal-safe colours
C_OVH = "#999999"  # grey
C_FLANK = "#6b9e6b"  # sage green
C_ARM = "#c9a94e"  # goldenrod
C_LOOP = "#b07ab0"  # muted purple
C_MM = "#ff6666"  # soft red — mismatch highlight
C_BB = "#5a7a9e"  # muted steel-blue — backbone dsDNA past the cut
C_SITE = "#4fa3a3"  # muted teal — recognition bases carried to rebuild the site

PAD = 3  # width of the "5' " strand-end markers
BB_EXT_LEN = 5  # how many "═" chars to extend the backbone past the cut
BB_CHAR = "═"  # dsDNA representation (double horizontal line)
GAP = 4  # blank columns between the backbone wings and the insert


def _split_highlighted(
    seq: str,
    default_colour: str,
    highlight_colour: str,
    positions: list[int],
) -> list[tuple[str, str]]:
    if not seq:
        return []
    pos_set = sorted({p for p in positions if 1 <= p <= len(seq)})
    if not pos_set:
        return [(default_colour, seq)]
    parts: list[tuple[str, str]] = []
    prev = 0
    for p in pos_set:
        idx = p - 1
        if idx > prev:
            parts.append((default_colour, seq[prev:idx]))
        parts.append((highlight_colour, seq[idx : idx + 1]))
        prev = idx + 1
    if prev < len(seq):
        parts.append((default_colour, seq[prev:]))
    return parts


def _bb_run() -> str:
    return f"[{C_BB}]{BB_CHAR * BB_EXT_LEN}[/{C_BB}]"


def _ss(seq: str) -> str:
    return f"[{C_OVH}]{seq}[/{C_OVH}]"


def _render_wing_left(cfg: OligoConfig) -> tuple[str, str, int]:
    """Backbone left arm, pre-ligation. Returns (top, bottom, width).

    Whichever strand of the *insert* protrudes at this junction, the
    backbone protrudes with the other one::

        5'-protruding insert  ->  backbone bottom sticks out right
        3'-protruding insert  ->  backbone top sticks out right

    A top line drawn left-to-right is the overhang's top-strand sequence
    verbatim; a bottom line is its positional complement.
    """
    ovh = cfg.five_prime_overhang
    bb = _bb_run()
    k = len(ovh)
    if ovh.kind == FIVE_PRIME:
        return f"{bb}{' ' * k}", f"{bb}{_ss(ovh.bottom)}", BB_EXT_LEN + k
    if ovh.kind == THREE_PRIME:
        return f"{bb}{_ss(ovh.top)}", f"{bb}{' ' * k}", BB_EXT_LEN + k
    return bb, bb, BB_EXT_LEN


def _render_wing_right(cfg: OligoConfig) -> tuple[str, str, int]:
    """Backbone right arm, pre-ligation. Returns (top, bottom, width).

    Mirror of `_render_wing_left`: a 5'-protruding insert end leaves the
    backbone's *top* strand sticking out to the left (e.g. XhoI's right arm
    really is ``5'-TCGAG...-3'``), and a 3'-protruding insert end leaves the
    backbone's bottom strand sticking out.
    """
    ovh = cfg.three_prime_overhang
    bb = _bb_run()
    k = len(ovh)
    if ovh.kind == FIVE_PRIME:
        return f"{_ss(ovh.top)}{bb}", f"{' ' * k}{bb}", k + BB_EXT_LEN
    if ovh.kind == THREE_PRIME:
        return f"{' ' * k}{bb}", f"{_ss(ovh.bottom)}{bb}", k + BB_EXT_LEN
    return bb, bb, BB_EXT_LEN


def core_parts(
    cfg: OligoConfig,
    arm1: str,
    arm2: str,
) -> list[tuple[str, str]]:
    """Colour-tagged segments of the insert's double-stranded core (5'->3')."""
    mm_positions = cfg.mismatch_positions

    # `mismatch_positions` are always counted from the base of the stem
    # (see `oligo.compute_mismatched_passenger`). arm1's front is always the
    # base of the stem (it opens onto the 5' flank) and arm2's front is
    # always the loop side (its back opens onto the 3' flank) — regardless
    # of which arm carries the guide vs. the passenger — so arm1 takes the
    # positions as-is and arm2 takes their mirror image about its own
    # length. Both arms are highlighted (not just the passenger's) since
    # the mismatch's paired base on the other strand forms the bulge too.
    arm1_mm = sorted({p for p in mm_positions if 1 <= p <= len(arm1)})
    arm2_mm = sorted(
        {len(arm2) - p + 1 for p in mm_positions if 1 <= len(arm2) - p + 1 <= len(arm2)}
    )

    # Mirrors `OligoConfig.core` segment for segment — the recognition bases
    # `keep_sites` restores bracket the flanks, and are empty otherwise.
    parts: list[tuple[str, str]] = [
        (C_SITE, cfg.site_bases(cfg.five_prime_overhang)),
        (C_FLANK, cfg.five_prime_flank),
    ]
    parts.extend(_split_highlighted(arm1, C_ARM, C_MM, arm1_mm))
    parts.append((C_LOOP, cfg.loop))
    parts.extend(_split_highlighted(arm2, C_ARM, C_MM, arm2_mm))
    parts.append((C_FLANK, cfg.three_prime_flank))
    parts.append((C_SITE, cfg.site_bases(cfg.three_prime_overhang)))
    return parts


def render_top_parts(
    cfg: OligoConfig, core: list[tuple[str, str]]
) -> list[tuple[str, str]]:
    """Top strand (5'->3'): core, plus whichever overhangs sit on the top.

    The left overhang joins the top strand when it is 5'-protruding, the
    right one when it is 3'-protruding.
    """
    left, right = cfg.five_prime_overhang, cfg.three_prime_overhang
    parts: list[tuple[str, str]] = []
    if left.kind == FIVE_PRIME:
        parts.append((C_OVH, left.top))
    parts.extend(core)
    if right.kind == THREE_PRIME:
        parts.append((C_OVH, right.top))
    return parts


def render_bottom_parts(
    cfg: OligoConfig, core: list[tuple[str, str]]
) -> list[tuple[str, str]]:
    """Bottom strand **as drawn left-to-right**, i.e. 3'->5'.

    Every column is the positional complement of the top-strand sequence
    there — including the columns where only the bottom strand exists,
    because overhangs are stored in top-strand coordinates. The left
    overhang joins the bottom strand when it is 3'-protruding, the right
    one when it is 5'-protruding.
    """
    left, right = cfg.five_prime_overhang, cfg.three_prime_overhang
    parts: list[tuple[str, str]] = []
    if left.kind == THREE_PRIME:
        parts.append((C_OVH, left.bottom))
    parts.extend((colour, complement(seq)) for colour, seq in core)
    if right.kind == FIVE_PRIME:
        parts.append((C_OVH, right.bottom))
    return parts


def render_alignment(
    top: str,
    top_offset: int,
    bottom: str,
    bottom_offset: int,
    width: int,
) -> str:
    """'|' wherever both strands are present and the bases pair."""
    out: list[str] = []
    for i in range(width):
        ti = i - top_offset
        bi = i - bottom_offset
        tb = top[ti] if 0 <= ti < len(top) else " "
        bb = bottom[bi] if 0 <= bi < len(bottom) else " "
        out.append("|" if COMPLEMENT.get(tb) == bb and bb != " " else " ")
    return "".join(out)


def _render_label_line(
    cfg: OligoConfig,
    arm1_len: int,
    arm2_len: int,
    width: int,
) -> str:
    """Segment-name annotations, in duplex column coordinates.

    The core always starts at column ``len(left overhang)``, so the two
    overhang labels bracket it regardless of which strand protrudes. Each
    overhang is labelled with its protrusion kind (``5'OVH`` / ``3'OVH``).
    """
    left, right = cfg.five_prime_overhang, cfg.three_prime_overhang

    def ovh_label(ovh: Overhang) -> str:
        # "5'OVH" wants 5 columns; fall back to the bare kind for the
        # 1- to 3-nt overhangs that NarI, BstNI and friends leave.
        return f"{ovh.kind[0]}'OVH" if len(ovh) >= 4 else f"{ovh.kind[0]}'"

    # Laid out with a running cursor rather than by summing the preceding
    # lengths at every entry: the segment list grew a `keep_sites` member at
    # each end, and hand-summed offsets are where that kind of change silently
    # slips a label one segment to the left.
    # The two `keep_sites` segments advance the cursor but carry no label: a
    # restriction-site remainder is usually a single base, and a name centred
    # on one column bleeds several columns either side, overwriting the
    # overhang and flank labels next to it. Their colour (`C_SITE`) is what
    # identifies them, and the caller states them in full alongside the
    # diagram.
    lengths = [
        (ovh_label(left), len(left)),
        ("", len(cfg.site_bases(left))),
        ("5'flank", len(cfg.five_prime_flank)),
        ("arm1", arm1_len),
        ("loop", len(cfg.loop)),
        ("arm2", arm2_len),
        ("3'flank", len(cfg.three_prime_flank)),
        ("", len(cfg.site_bases(right))),
        (ovh_label(right), len(right)),
    ]
    segments = []
    cursor = 0
    for name, length in lengths:
        segments.append((name, cursor, length))
        cursor += length

    buf_len = PAD + width + PAD
    chars: list[str] = [" "] * buf_len
    for name, start, length in segments:
        if length <= 0 or not name:
            continue
        centre = start + length // 2 + PAD
        l_start = max(0, centre - len(name) // 2)
        for j, ch in enumerate(name):
            pos = l_start + j
            if 0 <= pos < buf_len:
                chars[pos] = ch
    return "".join(chars).rstrip()


def _strand_row(
    rich: str,
    plain_len: int,
    offset: int,
    left_marker: str,
    right_marker: str,
    width: int,
) -> str:
    """Place a strand at ``offset`` in a ``PAD + width + PAD`` field.

    The end markers sit flush against the strand, so they stay meaningful
    when a strand starts or stops partway across the duplex — which is
    exactly what happens once one end is 3'-protruding.
    """
    pad = PAD
    lead = " " * (pad + offset - len(left_marker))
    trail = " " * (width - offset - plain_len + pad - len(right_marker))
    return f"{lead}{left_marker}{rich}{right_marker}{trail}"


def render_diagram(cfg: OligoConfig, arm1: str, arm2: str) -> list[str]:
    """The four-row junction diagram. Pure — no widget access."""
    core = core_parts(cfg, arm1, arm2)
    top_parts = render_top_parts(cfg, core)
    bot_parts = render_bottom_parts(cfg, core)
    duplex = make_duplex(cfg, arm1, arm2)
    width = duplex.width

    top_rich = "".join(f"[{c}]{s}[/{c}]" for c, s in top_parts if s)
    bot_rich = "".join(f"[{c}]{s}[/{c}]" for c, s in bot_parts if s)

    align = render_alignment(
        duplex.top, duplex.top_offset, duplex.bottom, duplex.bottom_offset, width
    )
    label_line = _render_label_line(cfg, len(arm1), len(arm2), width)

    left_top, left_bot, left_w = _render_wing_left(cfg)
    right_top, right_bot, right_w = _render_wing_right(cfg)

    pad, gap = PAD, GAP
    gap_sp = " " * gap
    left_pw = pad + left_w
    center_pw = pad + width + pad
    right_pw = right_w + pad
    total_w = left_pw + gap + center_pw + gap + right_pw

    # Rows 1/3: the two strands. The wings carry no end markers — the "═"
    # run already says the backbone continues past the frame.
    row1 = (
        f"{' ' * pad}{left_top}{gap_sp}"
        + _strand_row(top_rich, len(duplex.top), duplex.top_offset, "5' ", " 3'", width)
        + f"{gap_sp}{right_top}{' ' * pad}"
    )
    row2 = f"{' ' * left_pw}{gap_sp}{' ' * pad}{align}{' ' * pad}{gap_sp}{' ' * right_pw}"
    row3 = (
        f"{' ' * pad}{left_bot}{gap_sp}"
        + _strand_row(
            bot_rich, len(duplex.bottom), duplex.bottom_offset, "3' ", " 5'", width
        )
        + f"{gap_sp}{right_bot}{' ' * pad}"
    )

    # Row 4: wing captions plus the segment labels under the insert. The
    # buffer is oversized because a caption centred on a narrow wing runs
    # wider than the wing itself.
    buf_len = total_w + len("backbone right")
    chars = [" "] * buf_len
    for caption, centre in (
        ("backbone left", left_pw // 2),
        ("backbone right", left_pw + gap + center_pw + gap + right_pw // 2),
    ):
        for j, ch in enumerate(caption):
            pos = max(0, centre - len(caption) // 2) + j
            if pos < buf_len:
                chars[pos] = ch
    # Skip blanks so a caption that bleeds toward the insert is not erased
    # by the padding around the segment labels.
    for j, ch in enumerate(label_line):
        pos = left_pw + gap + j
        if ch != " " and pos < buf_len:
            chars[pos] = ch
    row4 = "".join(chars)

    return [row1.rstrip(), row2.rstrip(), row3.rstrip(), row4.rstrip()]
