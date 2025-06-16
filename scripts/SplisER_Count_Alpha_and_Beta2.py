import argparse
import gzip
import pandas as pd
import pybedtools
from pybedtools import BedTool
import sys
import psutil, os
import tempfile
import random
import string


# remove leftovers from previous crashes
pybedtools.cleanup(remove_all=True)

def print_memory_usage(msg=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024**2  # in MB
    print(f"[MEMORY] {msg} RSS: {mem:.2f} MB")

def bed12_juncs_to_bed6(bed12):
    for f in bed12:
        block_sizes = list(map(int, f.fields[10].split(',')))
        block_starts = list(map(int, f.fields[11].split(',')))
        if len(block_sizes) < 2 or len(block_starts) < 2:
            continue
        junc_start = f.start + block_starts[0] + block_sizes[0]
        junc_end = f.start + block_starts[1]
        yield pybedtools.create_interval_from_list([
            f.chrom, str(junc_start), str(junc_end), f.name, f.score, f.strand
        ])

def shorten_on_each_side(feature, len_to_shorten=1):
    """
    like bedtools slop -b but where b must be negative.
    """
    feature.start = feature.start + len_to_shorten
    feature.end = feature.end - len_to_shorten
    return feature

def try_delete_temporary_history(bedtool):
    try:
        bedtool.delete_temporary_history(ask=False)
    except AttributeError:
        pass

def reorder_bed_columns(df):
    # Reset index so 'name' is a column
    df = df.reset_index()
    # Define desired BED order (adjust as needed for your file)
    bed_cols = ['#Chrom', 'start', 'stop', 'name', 'TotalJuncCounts', 'strand']
    # Add sample columns (all columns not in bed_cols)
    sample_cols = [col for col in df.columns if col not in bed_cols]
    # Reorder
    df = df[bed_cols + sample_cols]
    return df

def count_alpha_beta2(splice_sites_df, junc_files, site_type):
    print(f"Counting {site_type} junctions for {len(junc_files)} files")
    alpha_counts = splice_sites_df.iloc[:,0:6].copy().set_index('name')
    beta2_counts = splice_sites_df.iloc[:,0:6].copy().set_index('name')

    # Make a set of splice sites for fast lookup
    site_set = set(splice_sites_df['name'].to_list())
    # Also, make splice sites as BedTool
    # Generate a random filename that does not match pybedtools.*.tmp
    rand_str = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    tmp_bed_path = f"/tmp/spliser_splice_sites_{rand_str}.bed"
    bedtool_splice_sites = pybedtools.BedTool(splice_sites_df.iloc[:,0:6].values.tolist()).sort().saveas(tmp_bed_path)


    # Prepare BedTool of all splice sites (as 1nt intervals)
    for junc_file in junc_files:
        print(f"Processing junction file: {junc_file}")
        sample = junc_file.split('/')[-1].replace('.junc', '')
        bed_juncs = BedTool(list(bed12_juncs_to_bed6(pybedtools.BedTool(junc_file))))
        # first calculate beta2 counts...
        beta2_counts_for_junc_file = dict()
        # make bed_juncs shortened by 1bp on each side so that we can use intersect to find encompassing junctions without including the junction itself
        bed_juncs_shortened = bed_juncs.each(shorten_on_each_side, len_to_shorten=1).sort()
        bed_splice_sites = bedtool_splice_sites.intersect(
            bed_juncs_shortened, s=True, wao=True, nonamecheck=True
        ).groupby(g=4, c=11, o='sum')
        # Read the groupby result into a dict
        beta2_counts_for_junc_file = {SpliceSite: 0 for SpliceSite in splice_sites_df['name'].to_list()}
        with open(bed_splice_sites.fn) as f:
            for line in f:
                SpliceSite, beta2_count = line.rstrip('\n').split('\t')
                if beta2_count == '-1':
                    continue
                else:
                    beta2_counts_for_junc_file[SpliceSite] = int(beta2_count)
        beta2_counts[sample] = pd.Series(beta2_counts_for_junc_file)

        # Now count alpha
        alpha_counts_for_junc_file={SpliceSite: 0 for SpliceSite in splice_sites_df['name'].to_list()}
        if site_type == 'Donors':
            for entry in bed_juncs:
                if entry.strand == '+':
                    SpliceSite_fromJunc = f'{entry.chrom}:{entry.start+1}:{entry.strand}:D'
                    if SpliceSite_fromJunc in site_set:
                        alpha_counts_for_junc_file[SpliceSite_fromJunc] += int(entry.score)
                else:
                    SpliceSite_fromJunc = f'{entry.chrom}:{entry.end}:{entry.strand}:D'
                    if SpliceSite_fromJunc in site_set:
                        alpha_counts_for_junc_file[SpliceSite_fromJunc] += int(entry.score)
        elif site_type == 'Acceptors':
            for entry in bed_juncs:
                if entry.strand == '+':
                    SpliceSite_fromJunc = f'{entry.chrom}:{entry.end}:{entry.strand}:A'
                    if SpliceSite_fromJunc in site_set:
                        alpha_counts_for_junc_file[SpliceSite_fromJunc] += int(entry.score)
                else:
                    SpliceSite_fromJunc = f'{entry.chrom}:{entry.start+1}:{entry.strand}:A'
                    if SpliceSite_fromJunc in site_set:
                        alpha_counts_for_junc_file[SpliceSite_fromJunc] += int(entry.score)
        alpha_counts[sample] = pd.Series(alpha_counts_for_junc_file)
        try_delete_temporary_history(bed_juncs_shortened)
        try_delete_temporary_history(bed_juncs)
        print_memory_usage(f"After processing {sample}")
        print(f"alpha_counts DataFrame size: {alpha_counts.memory_usage(deep=True).sum()/1024**2:.2f} MB")
        print(f"beta2_counts DataFrame size: {beta2_counts.memory_usage(deep=True).sum()/1024**2:.2f} MB")
        pybedtools.helpers.cleanup()
        # breakpoint()
        # print([(k,v) for (k,v) in alpha_counts_for_junc_file.items() if v > 0])
    alpha_counts = reorder_bed_columns(alpha_counts)
    beta2_counts = reorder_bed_columns(beta2_counts)
    return alpha_counts, beta2_counts

#i'm getting some occasional errors in
# bed_splice_sites = bedtool_splice_sites.intersect(bed_juncs_shortened, s=True, wao=True, nonamecheck =True).groupby(g=4, c=11, o='sum')

# and I think I would rather do the groupby sum just in pandas or however is most efficient in python, without using the pybedtools groupby. Like if we just take the bedtools_splice_sites.intersect(...) bedtool object, which is already sorted, by col4, and then sum based on col11 and use that to create beta2_counts_for_junc_file?

def main(argv=None):
    parser = argparse.ArgumentParser(description="Count Alpha and Beta2 junctions for each splice site.")
    parser.add_argument('--SpliceSites', required=True, help="Sorted splice site .bed.gz file")
    parser.add_argument('--InputJuncs', required=True, nargs='+', help="List of BED12 junc files")
    parser.add_argument('--AlphaOut', required=True, help="Output Alpha counts matrix (tsv.gz)")
    parser.add_argument('--Beta2Out', required=True, help="Output Beta2 counts matrix (tsv.gz)")
    parser.add_argument('--site_type', required=True, choices=["Donors", "Acceptors"], help="Splice site type to match")
    args = parser.parse_args(argv)
    splice_sites_df = pd.read_csv(args.SpliceSites, sep='\t')
    alpha_counts, beta2_counts = count_alpha_beta2(splice_sites_df, args.InputJuncs, args.site_type)

    alpha_counts.to_csv(args.AlphaOut, sep='\t', index=False)
    beta2_counts.to_csv(args.Beta2Out, sep='\t', index=False)

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        # short donors test
        main("--SpliceSites ../code/scratch/Donors.test.bed.gz --InputJuncs ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_1.junc ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_2.junc --AlphaOut ../code/scratch/alpha.tsv.gz --Beta2Out ../code/scratch/beta2.tsv.gz --site_type Donors".split())
        # short acceptors test
        main("--SpliceSites ../code/scratch/Acceptors.test.bed.gz --InputJuncs ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_1.junc ../code/rna-seq/SplicingAnalysis/juncfiles/X2_LCL_polyA_Branaplam_3.16_1.junc --AlphaOut ../code/scratch/alpha.tsv.gz --Beta2Out ../code/scratch/beta2.tsv.gz --site_type Acceptors".split())
        # longer donors test
        main("--SpliceSites ../code/rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.bed.gz --InputJuncs ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_1.junc ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_2.junc rna-seq/SplicingAnalysis/juncfiles/X4_LCL_polyA_WB01_NA_2.junc --AlphaOut ../code/scratch/alpha.tsv.gz --Beta2Out ../code/scratch/beta2.tsv.gz --site_type Donors".split())
        # longer acceptors test
        main("--SpliceSites ../code/rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Acceptors.bed.gz --InputJuncs ../code/rna-seq/SplicingAnalysis/juncfiles/1_Fibroblast_polyA_A10_NA_1.junc ../code/rna-seq/SplicingAnalysis/juncfiles/X4_LCL_polyA_WA11_NA_2.junc --AlphaOut ../code/scratch/alpha.tsv.gz --Beta2Out ../code/scratch/beta2.tsv.gz --site_type Acceptors".split())

    else:
        main()
