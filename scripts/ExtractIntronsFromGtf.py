import argparse
import gzip
import sys
import bedparse
import pysam
import subprocess
import tempfile
from Bio.Seq import Seq

def get_sequence(fasta, chrom, start, end, strand):
    seq = fasta.fetch(chrom, start, end)
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    return seq

def run_bedparse_gtf2bed(gtf_file, *args, n=None):
    """
    Wrapper around command line `bedparse gtf2bed`.
    Parameters:
    gtf_file (str): Path to the GTF file.
    *args: Additional arguments for the bedparse command.
    n (int, optional): Number of lines from the GTF file to process. If None, the entire file is processed.
    Returns:
    str: The stdout output from the bedparse command, or None if an error occurred.
    """
    # create temp file of first n lines
    if n is not None:
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_gtf:
            with open(gtf_file, 'r') as f:
                for i, line in enumerate(f):
                    if i >= n:
                        break
                    temp_gtf.write(line)
            gtf_file = temp_gtf.name
    # Construct the bedparse gtf2bed command
    command = ['bedparse', 'gtf2bed', gtf_file]
    command.extend(args)
    # Run the command and capture the output
    result = subprocess.run(command, capture_output=True, text=True)
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bedparse gtf2bed: {result.stderr}")
        return None
    return result.stdout

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Extract introns and their splice site sequences from GTF and reference FASTA.")
    parser.add_argument('--reference', required=True, help="Reference fasta file")
    parser.add_argument('--gtf', required=True, help="GTF file")
    parser.add_argument('--output', required=True, help="Output BED6 file")
    parser.add_argument('--n', type=int, help="Number of lines from the GTF file to process")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)
    
    # Extract introns from GTF file
    bed12 = run_bedparse_gtf2bed(args.gtf, n=args.n)
    if bed12 is None:
        return

    # Open the reference fasta file once
    fasta = pysam.FastaFile(args.reference)

    with open(args.output, 'w') as outfile:
        for i, l in enumerate(bed12.splitlines()):
            transcript = bedparse.bedline(l.split('\t')[0:12])
            introns = transcript.introns()
            if not introns:
                continue
            introns = introns.bed12tobed6()
            chrom = transcript.chr
            for intron in introns:
                start = intron.start
                end = intron.end
                strand = intron.strand
                # Extract 5'ss sequence (-4 to +7)
                if strand == '+':
                    five_ss_seq = get_sequence(fasta, chrom, start-4, start+7, strand)
                else:
                    five_ss_seq = get_sequence(fasta, chrom, end-7, end+4, strand)

                # Extract 3'ss sequence (-10 to +4)
                if strand == '+':
                    three_ss_seq = get_sequence(fasta, chrom, end-10, end+4, strand)
                else:
                    three_ss_seq = get_sequence(fasta, chrom, start-4, start+10, strand)

                # Write to output file in BED6 format
                name = f"{intron.name}|5'ss:{five_ss_seq}|3'ss:{three_ss_seq}"
                outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    # Close the reference fasta file
    fasta.close()

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        main("--reference /project2/yangili1/bjf79/ReferenceGenomes/Human_UCSC.hg38_GencodeComprehensive46/Reference.fa --gtf /project2/yangili1/bjf79/ReferenceGenomes/Human_UCSC.hg38_GencodeComprehensive46/Reference.gtf --output code/scratch/test.introns.bed --n 1000".split(' '))
    else:
        main()