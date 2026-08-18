#!/usr/bin/env bash
# Creates the small input files that .test/samples.tsv points at, so that a dry run finds every
# entry point present. Needs samtools on PATH.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data scratch

printf '@read1\nACGTACGTAC\n+\nIIIIIIIIII\n' | gzip > data/FromFastq.R1.fastq.gz
printf '@read1\nACGTACGTAC\n+\nIIIIIIIIII\n' | gzip > data/FromFastq.R2.fastq.gz

# chr1 length matches GRCh38, so this also passes the chromosome check against a real reference
{
    printf '@HD\tVN:1.6\tSO:coordinate\n'
    printf '@SQ\tSN:chr1\tLN:248956422\n'
    printf 'read1\t0\tchr1\t1000\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
} | samtools view -b -o data/FromBam.sorted.bam -

echo "wrote $(pwd)/data"
