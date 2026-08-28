#!/usr/bin/env python
"""
Merge per-sample compute-jui TSVs (rule ComputeJUI) into one junctions x samples wide matrix
(rule JUI_MergeMatrix), analogous to leafcutter's perind(_numers).counts.gz.

Engine: a streaming heapq.merge k-way merge across the already-sorted per-sample files. Chosen
over DuckDB/pandas/polars after benchmarking (see the project plan) because it holds at most one
junction's worth of rows in memory regardless of cohort size -- this rule reruns in full every
time a sample is added, so flat memory matters more than shaving off wall-clock time.
"""

import argparse
import ast
import csv
import gzip
import heapq
import itertools
import os

import pysam
from my_utils.jui_utils import DECOMPOSE_COLUMNS, TSV_COLUMNS

KNOWN_COLUMNS = set(TSV_COLUMNS) | set(DECOMPOSE_COLUMNS)
KEY_COLUMNS = ("chrom", "donor_pos", "acceptor_pos", "strand")

# Only these ast node types are permitted in --metric-expr -- blocks calls, attribute access,
# subscripting, comprehensions, etc. Keeps the expression to plain arithmetic over known columns.
_ALLOWED_NODES = (
    ast.Expression, ast.BinOp, ast.UnaryOp, ast.Name, ast.Load,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.USub, ast.UAdd,
    ast.Constant,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-tsvs", nargs="+", required=True, metavar="TSV",
                    help="Per-sample compute-jui output tsv.gz files (rule ComputeJUI). "
                         "Sample names are derived from each file's basename.")
    p.add_argument("--sample-suffix", default=".jui.tsv.gz",
                    help="Suffix stripped from each --input-tsvs basename to get the sample name "
                         "(default: .jui.tsv.gz, matching rule ComputeJUI's output naming).")
    p.add_argument("--metric-expr", required=True,
                    help="Arithmetic expression over compute-jui TSV columns, e.g. 'k' or 'R_D + R_A'.")
    p.add_argument("--fill", default="null",
                    help="Value for a sample that never observed a given junction: a number "
                         "(e.g. 0), or 'null'/'na'/'none' for NA (default: null). Note this only "
                         "governs true absence -- a junction a sample DID observe, where the "
                         "metric expression is itself undefined (e.g. a ratio below --min-count "
                         "in compute-jui), always renders as NA regardless of --fill.")
    p.add_argument("--min-total", type=float, default=None,
                    help="Drop junctions whose numeric values sum to less than this across all "
                         "samples (NA cells don't count toward the sum). Default: no filtering.")
    p.add_argument("--output", required=True, metavar="FILE.tsv.gz",
                    help="Output path. Written to FILE.tsv (uncompressed) first, then bgzipped "
                         "(and tabix-indexed if --tabix) to FILE.tsv.gz.")
    p.add_argument("--tabix", action="store_true",
                    help="bgzip + tabix-index --output in place, matching the per-sample files.")
    return p.parse_args()


def parse_fill(raw):
    if raw.strip().lower() in ("null", "na", "none", ""):
        return None
    return float(raw)


def compile_metric_expr(expr):
    tree = ast.parse(expr, mode="eval")
    names = set()
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise SystemExit(f"--metric-expr {expr!r}: disallowed syntax ({type(node).__name__}); "
                              f"only +, -, *, / and column names are permitted")
        if isinstance(node, ast.Name):
            names.add(node.id)
    unknown = names - KNOWN_COLUMNS
    if unknown:
        raise SystemExit(f"--metric-expr {expr!r}: unknown column(s) {sorted(unknown)}; "
                          f"known columns are {sorted(KNOWN_COLUMNS)}")
    return compile(tree, "<metric-expr>", "eval"), names


def sample_name_from_path(path, suffix):
    base = os.path.basename(path)
    if not base.endswith(suffix):
        raise SystemExit(f"{path}: basename does not end with --sample-suffix {suffix!r}")
    return base[: -len(suffix)]


def iter_sample_rows(path, sample_idx, needed_cols):
    """Yields (key, sample_idx, {col: raw_str}) for one sample's tsv.gz, key = (chrom, int(donor_pos), int(acceptor_pos), strand)."""
    with gzip.open(path, "rt", newline="") as fh:
        header = fh.readline()
        if not header.startswith("#"):
            raise SystemExit(f"{path}: expected a '#'-prefixed header line")
        cols = header[1:].rstrip("\n").split("\t")
        col_idx = {c: i for i, c in enumerate(cols)}
        missing = [c for c in (*KEY_COLUMNS, *needed_cols) if c not in col_idx]
        if missing:
            raise SystemExit(f"{path}: missing required column(s) {missing} (file has {cols})")
        ci = {c: col_idx[c] for c in (*KEY_COLUMNS, *needed_cols)}

        prev_key = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            key = (fields[ci["chrom"]], int(fields[ci["donor_pos"]]), int(fields[ci["acceptor_pos"]]), fields[ci["strand"]])
            if prev_key is not None and key < prev_key:
                raise SystemExit(
                    f"{path}: rows are not sorted by (chrom, donor_pos, acceptor_pos, strand) "
                    f"({prev_key} then {key}) -- MergeJUIMatrix requires compute-jui's native "
                    f"sorted output; do not reorder these files."
                )
            prev_key = key
            values = {c: fields[ci[c]] for c in needed_cols}
            yield key, sample_idx, values


def eval_metric(compiled_expr, values):
    """Returns a float, or None if any referenced column is NA/non-numeric for this row."""
    ns = {}
    for name, raw in values.items():
        if raw == "NA":
            return None
        try:
            ns[name] = float(raw)
        except ValueError:
            return None
    return eval(compiled_expr, {"__builtins__": {}}, ns)


def format_value(v):
    if v is None:
        return "NA"
    if v == int(v):
        return str(int(v))
    return str(round(v, 6))


_ABSENT = object()


def main():
    args = parse_args()
    fill = parse_fill(args.fill)
    compiled_expr, needed_cols = compile_metric_expr(args.metric_expr)

    samples = [sample_name_from_path(p, args.sample_suffix) for p in args.input_tsvs]
    if len(set(samples)) != len(samples):
        raise SystemExit(f"--input-tsvs basenames are not unique after stripping {args.sample_suffix!r}: {samples}")

    plain_path = args.output[:-3] if args.output.endswith(".gz") else args.output
    n = len(samples)

    generators = [
        iter_sample_rows(path, i, needed_cols)
        for i, path in enumerate(args.input_tsvs)
    ]
    merged = heapq.merge(*generators, key=lambda item: item[0])

    with open(plain_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["#chrom", "donor_pos", "acceptor_pos", "strand", "JunctionID", *samples])

        for key, group in itertools.groupby(merged, key=lambda item: item[0]):
            row = [_ABSENT] * n
            for _, sample_idx, values in group:
                row[sample_idx] = eval_metric(compiled_expr, values)

            resolved = [(fill if v is _ABSENT else v) for v in row]

            if args.min_total is not None:
                total = sum(v for v in resolved if v is not None)
                if total < args.min_total:
                    continue

            chrom, donor_pos, acceptor_pos, strand = key
            junction_id = f"{chrom}:{donor_pos}:{acceptor_pos}:{strand}"
            writer.writerow([chrom, donor_pos, acceptor_pos, strand, junction_id, *(format_value(v) for v in resolved)])

    if args.tabix:
        pysam.tabix_index(plain_path, seq_col=0, start_col=1, end_col=2, meta_char="#", zerobased=True, force=True)
    else:
        with open(plain_path, "rb") as src, gzip.open(args.output, "wb") as dst:
            dst.writelines(src)
        os.remove(plain_path)


if __name__ == "__main__":
    main()
