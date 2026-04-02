#!/usr/bin/env python3
"""
Comprehensive Sashimi-style plot for two samples with haplotypes.

Features:
- Two BAM/CRAM files (bam1, bam2)
- Split by HP tag (HP:1 / HP:2)
- Optional HP swap (--swap_hp bam1|bam2)
- Optional overlay of HP1+HP2 (--overlay bam1[,bam2])
- Optional highlight regions (--highlight "chr:start-end,chr:start-end")
- Optional highlight color (--highlight_color "#FFEE88" or "yellow")
- Optional GFF3 exon annotations (.gff3 or .gff3.gz)
- Custom sample labels (--l1, --l2)
- Figure size in pixels (--width, --height)
- Default width:height ratio = 2:1
- Resolution/DPI control (--res, default 300)
- Font and line thickness visually consistent; higher DPI = sharper output
- Scales Y-axis by 1.2 x clip_percentile of coverage (default 99)
"""

import argparse
import os
import gzip
from collections import Counter
import numpy as np
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl

# -----------------------
# Argument parsing
# -----------------------
def parse_args():
    p = argparse.ArgumentParser(description="Sashimi-style plot for two samples with haplotypes and optional GFF3.")
    p.add_argument("--bam1", required=True, help="Path to first BAM/CRAM file")
    p.add_argument("--bam2", required=True, help="Path to second BAM/CRAM file")
    p.add_argument("--region", required=True, help="Region to plot, format chr:start-end")
    p.add_argument("--highlight", default="", help="Comma-separated regions to highlight: chr:start-end,...")
    p.add_argument("--highlight_color", default="yellow",
                   help="Hex or named color for highlight (with or without #, default=yellow)")
    p.add_argument("--gff3", default=None, help="Optional GFF3 file (.gff3 or .gff3.gz)")
    p.add_argument("--swap_hp", default="", choices=["", "bam1", "bam2"],
                   help="Swap HP1/HP2 for one sample: 'bam1' or 'bam2'")
    p.add_argument("--overlay", default="", help="Comma-separated sample keys to overlay HPs in one panel: bam1,bam2")
    p.add_argument("--l1", default=None, help="Label for bam1 (default: bam1 filename)")
    p.add_argument("--l2", default=None, help="Label for bam2 (default: bam2 filename)")
    p.add_argument("--width", type=int, default=None, help="Figure width (pixels)")
    p.add_argument("--height", type=int, default=None, help="Figure height (pixels)")
    p.add_argument("--res", type=int, default=300, help="Output resolution (DPI), default 300")
    p.add_argument("--clip_percentile", type=float, default=99.0,
                   help="Percentile of coverage to use for Y-axis scaling (default 99)")
    p.add_argument("--output", default="sashimi_plot.png", help="Output image filename")
    return p.parse_args()

# -----------------------
# Helper functions
# -----------------------
def parse_gff3(gff3_file, chrom, start, end):
    exons = []
    if not gff3_file:
        return exons
    open_func = gzip.open if gff3_file.endswith(".gz") else open
    with open_func(gff3_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            c, _, typ, s, e, *_ = parts
            if c != chrom or typ != "exon":
                continue
            s, e = int(s), int(e)
            if e < start or s > end:
                continue
            exons.append((s, e))
    return exons

def parse_highlights(hstring):
    if not hstring:
        return []
    regions = []
    for seg in hstring.split(","):
        seg = seg.strip()
        if not seg:
            continue
        if ":" not in seg or "-" not in seg:
            raise ValueError(f"Invalid highlight region format: {seg}")
        regions.append(seg)
    return regions

def get_cov_junctions(bamfile, chrom, start, end, mapping_inv):
    sam = pysam.AlignmentFile(bamfile, "rb")
    cov_dict = {1:[0]*(end-start+1), 2:[0]*(end-start+1)}
    junctions_dict = {1:Counter(), 2:Counter()}

    for read in sam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("HP"):
            continue
        try:
            raw_hp = int(read.get_tag("HP"))
        except Exception:
            continue
        if raw_hp not in mapping_inv:
            continue
        panel_hp = mapping_inv[raw_hp]
        for pos in read.get_reference_positions():
            if start <= pos <= end:
                cov_dict[panel_hp][pos-start] += 1
        ref_pos = read.reference_start
        if read.cigartuples:
            for (ct, length) in read.cigartuples:
                if ct == 0:
                    ref_pos += length
                elif ct == 3:
                    junctions_dict[panel_hp][(ref_pos, ref_pos+length)] += 1
                    ref_pos += length
                elif ct == 2:
                    ref_pos += length
    sam.close()
    return cov_dict, junctions_dict

# -----------------------
# Plotting
# -----------------------
def plot_panel(ax, cov_dict, junctions_dict, start, end, colors,
               highlight_regions=None, exons=None, highlight_color='yellow', clip_percentile=99.0):
    if highlight_regions is None:
        highlight_regions = []
    if exons is None:
        exons = []

    x = range(start, end+1)
    max_height = 1

    for hp, cov in cov_dict.items():
        ax.fill_between(x, cov, step='mid', color=colors.get(hp,'gray'), alpha=0.5)
        if len(cov) > 0:
            # Clip by percentile to avoid extreme spikes
            if clip_percentile < 100:
                q = np.percentile(cov, clip_percentile)
            else:
                q = max(cov)
            max_height = max(max_height, q)

    for junctions in junctions_dict.values():
        for (jstart, jend), count in junctions.items():
            if jend < start or jstart > end:
                continue
            mid = (jstart+jend)/2
            ax.plot([jstart, mid, jend], [0, count, 0], color='red', alpha=0.7, linewidth=1)

    # --- Highlights ---
    hl_color = highlight_color.strip()
    def looks_like_hex(s):
        return all(c in "0123456789ABCDEFabcdef" for c in s)
    if not hl_color.startswith("#") and looks_like_hex(hl_color) and len(hl_color) in (3, 6):
        hl_color = f"#{hl_color}"

    for r in highlight_regions:
        chr_, pos = r.split(":")
        hs, he = map(int, pos.split("-"))
        if he < start or hs > end:
            continue
        ax.axvspan(max(hs, start), min(he, end), color=hl_color, alpha=0.3)

    # --- Exons ---
    exon_y_bottom = max_height * 0.02
    exon_height = max(1.0, max_height * 0.05)
    for s, e in exons:
        if e < start or s > end:
            continue
        s, e = max(s, start), min(e, end)
        ax.fill_between([s, e], [exon_y_bottom, exon_y_bottom],
                        [exon_y_bottom+exon_height, exon_y_bottom+exon_height],
                        color='green', alpha=0.6, linewidth=0)

    ax.set_xlim(start, end)
    ax.set_ylim(0, max_height * 1.2)
    ax.set_ylabel("Coverage")
    ax.grid(True, linestyle='--', alpha=0.3)

# -----------------------
# Main
# -----------------------
def main():
    args = parse_args()

    # Parse region
    if ":" not in args.region or "-" not in args.region:
        raise ValueError("Region must be chr:start-end")
    chrom, pos = args.region.split(":")
    start, end = map(int, pos.split("-"))
    if end <= start:
        raise ValueError("End must be greater than start")

    # Highlights
    highlights = parse_highlights(args.highlight)
    highlight_color = args.highlight_color.strip()

    # Overlay samples
    overlay_set = set([s.strip() for s in args.overlay.split(",") if s.strip()])

    # Labels
    label1 = args.l1 if args.l1 else os.path.basename(args.bam1)
    label2 = args.l2 if args.l2 else os.path.basename(args.bam2)
    samples = [("bam1", label1, args.bam1), ("bam2", label2, args.bam2)]

    # Determine panel count
    n_panels = sum(1 if key in overlay_set else 2 for key, _, _ in samples)

    # --- Figure size and DPI ---
    default_w, default_h = 1600, 800
    width_px = args.width or (args.height * 2 if args.height else default_w)
    height_px = args.height or (args.width / 2 if args.width else default_h)
    dpi = args.res if args.res > 0 else 300

    fig_w_in = width_px / 100.0
    fig_h_in = height_px / 100.0

    # --- Keep fonts & lines visually consistent ---
    rcp = mpl.rcParams
    rcp['figure.dpi'] = 100
    rcp['savefig.dpi'] = dpi

    # Exons
    exons = parse_gff3(args.gff3, chrom, start, end) if args.gff3 else []
    colors = {1:'skyblue', 2:'lightcoral'}

    # Create figure
    fig, axes = plt.subplots(n_panels, 1, figsize=(fig_w_in, fig_h_in), sharex=True)
    if n_panels == 1:
        axes = [axes]

    panel_idx = 0
    for key, label, bam in samples:
        mapping_inv = {1:2,2:1} if args.swap_hp == key else {1:1,2:2}
        cov_dict, junc_dict = get_cov_junctions(bam, chrom, start, end, mapping_inv)
        if key in overlay_set:
            ax = axes[panel_idx]
            plot_panel(ax, cov_dict, junc_dict, start, end, colors,
                       highlight_regions=highlights, exons=exons,
                       highlight_color=highlight_color, clip_percentile=args.clip_percentile)
            ax.set_title(f"{label} (HP1 + HP2)")
            panel_idx += 1
        else:
            for hp in (1,2):
                ax = axes[panel_idx]
                plot_panel(ax, {hp:cov_dict[hp]}, {hp:junc_dict[hp]}, start, end, colors,
                           highlight_regions=highlights, exons=exons,
                           highlight_color=highlight_color, clip_percentile=args.clip_percentile)
                ax.set_title(f"{label} HP{hp}")
                panel_idx += 1

    axes[-1].set_xlabel(f"{chrom}:{start}-{end}")
    plt.tight_layout()
    plt.savefig(args.output)
    plt.close(fig)
    print(f" Sashimi plot saved to {args.output} ({width_px}x{height_px}px at {dpi} DPI)")
    print(f"   Y-scale based on {args.clip_percentile}th percentile of coverage per haplotype")

if __name__ == "__main__":
    main()

