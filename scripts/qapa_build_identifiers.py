#!/usr/bin/env python3
"""
Parse a Gencode/Ensembl GTF to build the tab-delimited identifiers table
expected by qapa (--db flag). Reads transcript records only, deduplicated.

Usage: python qapa_build_identifiers.py input.gtf output.txt
"""
import sys
import re


def parse_attr(attr_str, key):
    m = re.search(rf'{key} "([^"]+)"', attr_str)
    return m.group(1) if m else ""


def strip_version(ensembl_id):
    """Strip Ensembl version suffix: ENST00000456328.2 -> ENST00000456328.
    QAPA's get_stripped_name() does this before looking up the identifiers table,
    so IDs in the table must not have version suffixes."""
    return ensembl_id.split(".")[0]


def main():
    gtf_file = sys.argv[1]
    out_file = sys.argv[2]

    seen = set()
    with open(gtf_file) as fh, open(out_file, "w") as out:
        out.write(
            "Gene stable ID\tTranscript stable ID\t"
            "Gene type\tTranscript type\tGene name\n"
        )
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            gene_id = strip_version(parse_attr(attrs, "gene_id"))
            transcript_id = strip_version(parse_attr(attrs, "transcript_id"))
            gene_type = parse_attr(attrs, "gene_type") or parse_attr(
                attrs, "gene_biotype"
            )
            transcript_type = parse_attr(
                attrs, "transcript_type"
            ) or parse_attr(attrs, "transcript_biotype")
            gene_name = parse_attr(attrs, "gene_name")
            key = (gene_id, transcript_id)
            if key in seen:
                continue
            seen.add(key)
            out.write(
                f"{gene_id}\t{transcript_id}\t{gene_type}\t"
                f"{transcript_type}\t{gene_name}\n"
            )


if __name__ == "__main__":
    main()
