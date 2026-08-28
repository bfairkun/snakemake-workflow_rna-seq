#!/usr/bin/env python
"""
Scale a bigwig to a target total covered bases genome-wide, reading and writing via pyBigWig.

Replaces the old idxstats-read-count-based normalization (bigWigToBedGraph | awk | sort |
bedGraphToBigWig): the scale factor's denominator is now total bases covered (sum of
depth * interval-width), read straight from the bigwig's own data rather than a separate
samtools idxstats pass over the bam. Total bases covered is the more principled denominator
for scaling a coverage track, since a bigwig's units are exactly depth * width.

Processes one fixed-size genomic window at a time (--chunk-size), so memory stays bounded by
window size and interval density regardless of chromosome/genome size -- the whole normalized
bigwig is never held in memory at once.
"""

import argparse
import re
from math import floor, log10

import pyBigWig


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, metavar="RAW.bw", help="Raw (unnormalized) input bigwig.")
    p.add_argument("--output", required=True, metavar="NORMALIZED.bw", help="Normalized output bigwig.")
    p.add_argument(
        "--chrom-filter", nargs="+", default=["^chr[0-9]+$"], metavar="REGEX",
        help="Regex(es) (any match counts) selecting which chromosomes' covered bases sum to "
             "the scale-factor denominator. Does NOT restrict which chromosomes are written to "
             "the output -- every chromosome in the input bigwig is written, all scaled by the "
             "same single factor. Default: numbered chromosomes only (excludes chrX/Y/M and any "
             "non-'chrN' contigs).",
    )
    p.add_argument(
        "--target-total-bases", type=float, default=1_000_000_000, metavar="N",
        help="Scale so that total covered bases (on --chrom-filter-matching chromosomes) equals "
             "this value (default: 1e9, i.e. coverage per billion covered bases).",
    )
    p.add_argument(
        "--chunk-size", type=int, default=10_000_000, metavar="BP",
        help="Genomic window size in bp processed per read/scale/write step (default: 10,000,000). "
             "Bounds memory use independent of chromosome/genome size; does not affect output values.",
    )
    p.add_argument(
        "--sig-digits", type=int, default=5, metavar="N",
        help="Round each scaled value to this many significant figures before writing (default: 5). "
             "Keeps file size down; the original awk-based pipeline had ~6 via awk's default OFMT. "
             "Set to 0 to disable rounding.",
    )
    return p.parse_args()


def round_sig(x, sig):
    if sig <= 0 or x == 0:
        return x
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def total_covered_bases(bw, chroms, patterns):
    matching = [c for c in chroms if any(p.search(c) for p in patterns)]
    if not matching:
        raise SystemExit(
            f"No chromosomes in the input bigwig matched any of {[p.pattern for p in patterns]} "
            f"(bigwig has: {sorted(chroms)}); scale factor would be undefined."
        )
    total = 0.0
    for c in matching:
        s = bw.stats(c, 0, chroms[c], type="sum", exact=True)[0]
        total += s or 0.0
    if total <= 0:
        raise SystemExit(
            f"Total covered bases on chromosomes matching {[p.pattern for p in patterns]} is zero; "
            f"scale factor would be undefined."
        )
    return total


def main():
    args = parse_args()
    patterns = [re.compile(p) for p in args.chrom_filter]

    bw_in = pyBigWig.open(args.input)
    chroms = bw_in.chroms()  # {name: length}, also gives us chrom sizes for the header -- no .fai needed

    total_bases = total_covered_bases(bw_in, chroms, patterns)
    scale_factor = args.target_total_bases / total_bases

    bw_out = pyBigWig.open(args.output, "w")
    bw_out.addHeader(list(chroms.items()))

    for chrom, length in chroms.items():
        for window_start in range(0, length, args.chunk_size):
            window_end = min(window_start + args.chunk_size, length)
            intervals = bw_in.intervals(chrom, window_start, window_end)
            if not intervals:
                continue
            # bw.intervals() returns intervals overlapping the query window unclipped, so an
            # interval straddling a chunk boundary would otherwise be re-emitted (unclipped) in
            # both neighboring windows -- clip to this window to keep chunks non-overlapping.
            starts = [max(i[0], window_start) for i in intervals]
            ends = [min(i[1], window_end) for i in intervals]
            values = [round_sig(i[2] * scale_factor, args.sig_digits) for i in intervals]
            bw_out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw_out.close()
    bw_in.close()


if __name__ == "__main__":
    main()
