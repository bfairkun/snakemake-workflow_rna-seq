import argparse
import pandas as pd
import pysam
from Bio import motifs
from Bio.Seq import Seq
import sys

def get_sequence(fasta, chrom, start, end, strand):
    seq = fasta.fetch(chrom, start, end)
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    return seq

def parse_introns_with_ss(introns_file):
    introns = pd.read_csv(introns_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'ss', 'score', 'strand'])
    introns['5ss'] = introns['ss'].apply(lambda x: x.split('|')[0].split(':')[1])
    introns['3ss'] = introns['ss'].apply(lambda x: x.split('|')[1].split(':')[1])
    return introns

def create_pwm(sequences):
    instances = [Seq(seq) for seq in sequences]
    m = motifs.create(instances)
    pwm = m.counts.normalize(pseudocounts=0.5)
    return pwm

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Add splice site sequences and scores to regtools annotate output.")
    parser.add_argument('--input', required=True, help="Input annotated junctions file")
    parser.add_argument('--reference', required=True, help="Reference fasta file")
    parser.add_argument('--introns', required=True, help="Annotated introns with splice site sequences file")
    parser.add_argument('--output', required=True, help="Output file")
    parser.add_argument('--nrows', type=int, default=None, help="Number of rows to read from input file (optional)")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    # Read input files
    annotated_junctions = pd.read_csv(args.input, sep='\t', nrows=args.nrows)
    introns_with_ss = parse_introns_with_ss(args.introns)

    # Create PWMs for 5'ss and 3'ss
    pssm_5ss = create_pwm(introns_with_ss['5ss']).log_odds()
    pssm_3ss = create_pwm(introns_with_ss['3ss']).log_odds()

    # Open the reference fasta file once
    fasta = pysam.FastaFile(args.reference)

    # annotated_junctions = annotated_junctions.head(1000) # For testing

    # Add new columns to annotated_junctions
    annotated_junctions['5ss_seq'] = annotated_junctions.apply(
        lambda row: get_sequence(fasta, row['chrom'], row['end']-8, row['end']+3, row['strand']) if row['strand'] == '-' else get_sequence(fasta, row['chrom'], row['start']-4, row['start']+7, row['strand']), axis=1)
    annotated_junctions['3ss_seq'] = annotated_junctions.apply(
        lambda row: get_sequence(fasta, row['chrom'], row['start']-4, row['start']+10, row['strand']) if row['strand'] == '-' else get_sequence(fasta, row['chrom'], row['end']-11, row['end']+3, row['strand']), axis=1)

    annotated_junctions['5ss_score'] = annotated_junctions['5ss_seq'].apply(lambda seq: pssm_5ss.calculate(seq))
    annotated_junctions['3ss_score'] = annotated_junctions['3ss_seq'].apply(lambda seq: pssm_3ss.calculate(seq))

    # Write output file
    annotated_junctions.to_csv(args.output, sep='\t', index=False)

    # Close the reference fasta file
    fasta.close()

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        main("--input rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/GRCh38_GencodeRelease44Comprehensive.uniq.annotated.tsv.gz --reference /project2/yangili1/bjf79/ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Reference.fa --introns /project2/yangili1/bjf79/ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Reference.Introns.bed.gz --output scratch/GenomeName.uniq.annotated.with_ss_scores.tsv.gz --nrows 10000".split(' '))
    else:
        main()