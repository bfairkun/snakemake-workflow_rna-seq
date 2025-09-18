import argparse
import pandas as pd
import numpy as np
import rnanorm
import re

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert featureCounts output to log2(TPM), log2(TMM-CPM), and log2(filtered TMM-CPM) matrices in BED6+ format."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input featureCounts count file (tab-separated, with '#' comment line)"
    )
    parser.add_argument(
        "--log2TPM_Matrix_bed", required=True, help="Output BED6+ file for log2(TPM) matrix"
    )
    parser.add_argument(
        "--log2TMM_Normalized_CPM_Matrix_bed", required=True, help="Output BED6+ file for log2(TMM-normalized CPM) matrix"
    )
    parser.add_argument(
        "--log2Filtered_TMM_Normalized_CPM_Matrix_bed", required=True, help="Output BED6+ file for log2(filtered TMM-normalized CPM) matrix"
    )
    parser.add_argument(
        "--pseudocount", type=float, default=1.0, help="Pseudocount to add before log transform (default: 1)"
    )
    parser.add_argument(
        "--cpm_filter", type=float, default=1.0, help="Mean CPM threshold for filtering genes (default: 1)"
    )
    parser.add_argument(
        "--sample_rename_regex",
        type=str,
        default=None,
        help="Optional regex with one capture group to rename sample columns, e.g. 'Alignments/(.*?)/Aligned'"
    )
    return parser.parse_args(argv)

def read_featurecounts(file):
    # Find header line (starts with Geneid)
    with open(file) as f:
        for i, line in enumerate(f):
            if (line.startswith("Geneid")):
                header = line.strip().split('\t')
                skip = i
                break
    df = pd.read_csv(file, sep='\t', comment='#', skiprows=skip, header=0)
    return df

def get_bed6(df):
    def collapse_chrom(x):
        return str(x.split(';')[0]) if ';' in x else x
    def collapse_start(x):
        return int(min([int(i) for i in x.split(';')]))
    def collapse_end(x):
        return int(max([int(i) for i in x.split(';')]))
    def collapse_strand(x):
        return x.split(';')[0] if ';' in x else x

    bed6 = pd.DataFrame({
        '#chrom': df['Chr'].apply(collapse_chrom),
        'start': df['Start'].apply(collapse_start).astype(int),
        'end': df['End'].apply(collapse_end).astype(int),
        'name': df['Geneid'] if 'Geneid' in df.columns else df.index,
        'score': 0,
        'strand': df['Strand'].apply(collapse_strand),
    })
    # Ensure start and end are int dtype
    bed6['start'] = bed6['start'].astype(int)
    bed6['end'] = bed6['end'].astype(int)
    return bed6

def main(argv=None):
    args = parse_args(argv)
    df = read_featurecounts(args.input).set_index("Geneid")
    sample_cols = [col for col in df.columns if col not in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']]

    # Optional: Rename sample columns using regex
    if args.sample_rename_regex:
        regex = re.compile(args.sample_rename_regex)
        new_sample_cols = []
        for col in sample_cols:
            m = regex.fullmatch(col)  # <-- use fullmatch instead of search
            if m:
                new_sample_cols.append(m.group(1))
            else:
                new_sample_cols.append(col)
        # Rename columns in counts and all downstream DataFrames
        df = df.rename(columns=dict(zip(sample_cols, new_sample_cols)))
        sample_cols = new_sample_cols

    # Prepare counts and gene lengths
    counts = df[sample_cols].astype(float)  # Geneid as index
    gene_lengths = df['Length'].astype(float)

    # TPM normalization
    tpm = rnanorm.TPM(gene_lengths=gene_lengths).fit_transform(counts.T + args.pseudocount)
    tpm = pd.DataFrame(tpm, index=counts.columns, columns=counts.index).T  # genes x samples
    log2tpm = np.log2(tpm)
    bed6 = get_bed6(df)
    bed6['score'] = log2tpm.mean(axis=1).values
    # Ensure indices match for concat
    log2tpm_out = log2tpm[sample_cols]
    log2tpm_out.index = bed6.index  # align indices
    log2tpm_df = pd.concat([bed6, log2tpm_out], axis=1)
    log2tpm_df.to_csv(args.log2TPM_Matrix_bed, sep='\t', index=False, float_format='%.4f')

    # CPM calculation
    lib_sizes = counts.sum(axis=0)
    cpm = (counts / lib_sizes) * 1e6

    # TMM normalization on CPM
    tmm_cpm = rnanorm.TMM().fit_transform(counts.T + args.pseudocount)
    tmm_cpm = pd.DataFrame(tmm_cpm, index=counts.columns, columns=counts.index).T
    log2tmm_cpm = np.log2(tmm_cpm)
    bed6['score'] = log2tmm_cpm.mean(axis=1).values
    # Ensure indices match for concat
    log2tmm_cpm_out = log2tmm_cpm[sample_cols]
    log2tmm_cpm_out.index = bed6.index  # align indices
    log2tmm_cpm_df = pd.concat([bed6, log2tmm_cpm_out], axis=1)
    log2tmm_cpm_df.to_csv(args.log2TMM_Normalized_CPM_Matrix_bed, sep='\t', index=False, float_format='%.4f')

    # Filter genes with mean CPM > threshold
    mean_cpm = cpm.mean(axis=1)
    filter_mask = mean_cpm > args.cpm_filter
    filtered_counts = counts[filter_mask]
    filtered_df = df[filter_mask].reset_index()  # Do NOT use drop=True
    filtered_bed6 = get_bed6(filtered_df)

    # TMM normalization on filtered CPM
    filtered_tmm_cpm = rnanorm.TMM().fit_transform(filtered_counts.T + args.pseudocount)
    filtered_tmm_cpm = pd.DataFrame(filtered_tmm_cpm, index=filtered_counts.columns, columns=filtered_counts.index).T
    log2filtered_tmm_cpm = np.log2(filtered_tmm_cpm)
    filtered_bed6['score'] = log2filtered_tmm_cpm.mean(axis=1).values
    # Ensure indices match for concat
    log2filtered_tmm_cpm_out = log2filtered_tmm_cpm[sample_cols]
    log2filtered_tmm_cpm_out.index = filtered_bed6.index  # align indices
    log2filtered_tmm_cpm_df = pd.concat([filtered_bed6, log2filtered_tmm_cpm_out], axis=1)
    log2filtered_tmm_cpm_df.to_csv(args.log2Filtered_TMM_Normalized_CPM_Matrix_bed, sep='\t', index=False, float_format='%.4f')

if __name__ == "__main__":
    import sys
    if hasattr(sys, 'ps1') or sys.flags.interactive:
        # Example: main("your arguments here".split())
        main("-i rna_seq/featureCounts/GRCh38_GencodeRelease44Comprehensive/AllSamplesUnstrandedCounting.Counts.txt --log2TPM_Matrix_bed scratch/log2TPM.bed --log2TMM_Normalized_CPM_Matrix_bed scratch/log2TMM_CPM.bed --log2Filtered_TMM_Normalized_CPM_Matrix_bed scratch/log2Filtered_TMM_CPM.bed".split())
    else:
        main()