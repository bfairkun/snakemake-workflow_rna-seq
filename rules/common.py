import pandas as pd
import numpy as np
import os
import shutil
import subprocess
import requests
import yaml

from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.io import glob_wildcards
from snakemake.exceptions import WorkflowError



min_version("7.32")

# Define default configuration for interactive use
if 'config' not in globals():
    config = {
        "samples": "config/samples.tsv",
        "STAR_genomes": "module_workflows/rna_seq/config/STAR_Genome_List.tsv",
        "contrast_group_files_prefix": "config/contrast_group_files/"
    }

# general functions
def read_config_file(config_file):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def remove_prefix(s, prefix):
    return s[len(prefix):] if s.startswith(prefix) else s

def fetch_ena_links(accession):
    """Fetch R1 and R2 FTP links for a given SRA accession."""
    if (accession in [np.nan, ""]) or pd.isna(accession):
        # print("None!!")
        return None, None
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "json"
    }
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        if not data:
            return None, None  # No links found
        fastq_links = data[0]["fastq_ftp"]
        links = ["ftp://" + link for link in fastq_links.split(";")]  # Split R1 and R2 links
        return links[0], links[1] if len(links) > 1 else None
    except Exception as e:
        print(f"Error fetching links for {accession}: {e}")
        return None, None

def is_given(value):
    return pd.notna(value) and value != ""

def populate_Library_Layout(row):
    if is_given(row['R1']) and is_given(row['R2']):
        return "PAIRED"
    elif is_given(row['R1']):
        return "SINGLE"
    elif is_given(row['R1_link']) and is_given(row['R2_link']):
        return "PAIRED"
    elif is_given(row['R1_link']):
        return "SINGLE"
    else:
        return ""

def read_fai_lengths(fai):
    lengths = {}
    with open(fai) as f:
        for line in f:
            fields = line.split('\t')
            lengths[fields[0]] = int(fields[1])
    return lengths

def read_bam_header(bam):
    """Returns (is_coordinate_sorted, {chrom: length}) from the bam header"""
    header = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
    is_sorted = False
    lengths = {}
    for line in header.split('\n'):
        fields = line.split('\t')
        if fields[0] == '@HD':
            is_sorted = 'SO:coordinate' in fields
        elif fields[0] == '@SQ':
            tags = dict(tag.split(':', 1) for tag in fields[1:] if ':' in tag)
            if 'SN' in tags and 'LN' in tags:
                lengths[tags['SN']] = int(tags['LN'])
    return is_sorted, lengths

def compare_chroms_to_reference(bam_lengths, fai_lengths):
    """Returns None if the bam is compatible with the reference, else a reason string"""
    shared = set(bam_lengths) & set(fai_lengths)
    if not shared:
        return (f"no chromosome names in common with the reference "
                f"(bam has {sorted(bam_lengths)[:3]}, reference has {sorted(fai_lengths)[:3]}); "
                f"check for a chr-prefix naming mismatch")
    mismatched = [c for c in sorted(shared) if bam_lengths[c] != fai_lengths[c]]
    if mismatched:
        c = mismatched[0]
        return (f"chromosome {c} is {bam_lengths[c]} bp in the bam but {fai_lengths[c]} bp in the "
                f"reference ({len(mismatched)} mismatched in total); likely a different assembly")
    return None

def precheck_bam(bam, fai):
    """Returns None if the bam is usable, else a reason string. Chromosome names are only
    checked if the reference has already been indexed; rule StageProvidedBam always checks."""
    if not os.path.exists(bam):
        return "file not found"
    try:
        is_sorted, bam_lengths = read_bam_header(bam)
    except subprocess.CalledProcessError as e:
        return f"samtools could not read the header: {e.stderr.strip()}"
    if not is_sorted:
        return "header lacks @HD SO:coordinate; bam must be coordinate-sorted"
    if os.path.exists(fai):
        return compare_chroms_to_reference(bam_lengths, read_fai_lengths(fai))
    return None

def wildcard_constraints_from_list(l):
    if len(l) == 0:
        return "DUMMY"
    else:
        return "|".join(l)

# read and process samples files
samples = pd.read_csv(config["samples"],sep='\t')
STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0, converters={"ChromLargerThan512Mbp": lambda x: x.lower() in ["yes", "true", "y", "True", "T"]})

STAR_genomes['GenomeName'] = STAR_genomes.index

# Add optional columns if not existing
for optional_column in ['bam', 'R1', 'R2', 'R1_link', 'R2_link', 'SRA_accession', 'Library_Layout']:
    if optional_column not in samples.columns:
        samples[optional_column] = pd.NA

# A provided bam takes precedence over R1/R2 and SRA_accession, but only if it is usable.
# If it isn't, fall back to fastq for that sample when a fastq source is also provided, since
# the entry point has to be settled here, before the DAG is built.
if config.get('bam_precheck', True) and samples['bam'].apply(is_given).any():
    if shutil.which("samtools") is None:
        print("Warning: samtools not found on PATH, skipping bam prechecks. Provided bams will still be checked by rule StageProvidedBam.")
    else:
        unusable = []
        for i, row in samples.loc[samples['bam'].apply(is_given)].iterrows():
            reason = precheck_bam(row['bam'], config['GenomesPrefix'] + row['STARGenomeName'] + "/Reference.fa.fai")
            if reason is None:
                continue
            if is_given(row['R1']) or is_given(row['R1_link']) or is_given(row['SRA_accession']):
                print(f"Warning: sample {row['sample']} bam {row['bam']}: {reason}. Falling back to the fastq/SRA source given in the same row.")
                samples.loc[i, 'bam'] = pd.NA
            else:
                unusable.append(f"  {row['sample']}: {row['bam']}\n    {reason}")
        if unusable:
            raise WorkflowError(
                "Unusable bam file(s) with no fastq/SRA fallback in the same row:\n" + "\n".join(unusable) +
                "\nFix the bam(s), or provide R1/R2 or SRA_accession in the same row to fall back to aligning from fastq."
            )

samples['bam'] = samples['bam'].where(samples['bam'].apply(is_given), pd.NA)

# Fill fastq links from SRA accession, except where a usable bam makes them unnecessary
samples[["R1_link", "R2_link"]] = samples.apply(
    lambda row: pd.Series(fetch_ena_links(row["SRA_accession"]))
    if pd.isna(row["bam"]) and (pd.isna(row["R1_link"]) or row["R1_link"] == "") and not pd.isna(row["SRA_accession"])
    else pd.Series([row["R1_link"], row["R2_link"]]),
    axis=1
)

# Add column for SE vs PE if contains any NAs/blanks or column doesn't already exist
samples['Library_Layout'] = samples.apply(lambda row: populate_Library_Layout(row) if pd.isna(row['Library_Layout']) or row['Library_Layout'] == "" else row['Library_Layout'], axis=1)

samples_missing_layout = samples.loc[samples['bam'].notna() & ~samples['Library_Layout'].apply(is_given)]['sample'].unique()
if len(samples_missing_layout) > 0:
    raise WorkflowError(
        f"Library_Layout (PAIRED or SINGLE) must be given for samples entering the workflow as a bam, "
        f"since it cannot be inferred from R1/R2: {', '.join(samples_missing_layout)}"
    )

# validate() only inserts a schema default for a column that is absent altogether, so blank cells
# in a column that does exist keep their NaN and end up in filenames as "nan". Fill them first.
for column, spec in read_config_file(os.path.join(workflow.current_basedir, "../schemas/samples.schema.yaml"))['properties'].items():
    if 'default' in spec and column in samples.columns:
        samples[column] = samples[column].where(samples[column].apply(is_given), spec['default'])

# validate
validate(samples, "../schemas/samples.schema.yaml")
validate(STAR_genomes, "../schemas/STAR_Genome_List.schema.yaml")

AllSamples = samples['sample'].unique()

samples_SingleEnd = samples[samples['Library_Layout']=='SINGLE']['sample']
samples_PairedEnd = samples[samples['Library_Layout']=='PAIRED']['sample']

samples_from_links = samples.loc[samples['bam'].isna() & samples["R1"].isna() & samples["R1_link"].notna()]['sample']
samples_from_local = samples.loc[samples['bam'].isna() & samples["R1"].notna()]['sample']
samples_from_bam = samples.loc[samples['bam'].notna()]['sample'].unique()

samples_mixing_bam_and_fastq = [s for s in samples_from_bam if s in set(samples_from_local) | set(samples_from_links)]
if samples_mixing_bam_and_fastq:
    raise WorkflowError(
        f"A sample must enter the workflow either as bam(s) or as fastq, not both. Rows for these samples "
        f"provide a bam and fastq/SRA sources for the same sample: {', '.join(samples_mixing_bam_and_fastq)}"
    )

# Samples entering as a bam are already aligned; Aligner still describes what kind of reads they
# are, so they remain included in the short-read (STAR) quantification rules
samples_needing_alignment = [s for s in AllSamples if s not in set(samples_from_bam)]

samples_ForSTAR_df = samples.loc[samples['Aligner']=='STAR']
samples_ForSTAR = samples_ForSTAR_df['sample'].unique()
samples_STAR_aligned_here = samples_ForSTAR_df.loc[samples_ForSTAR_df['bam'].isna()]['sample'].unique()

SingleEndSamples_wildcards_regex = wildcard_constraints_from_list(samples_SingleEnd)
PairedEndSamples_wildcards_regex = wildcard_constraints_from_list(samples_PairedEnd)

# differential splicing and expression contrasts
if config["contrast_group_files_prefix"]:
    contrasts, = glob_wildcards(f"{config['contrast_group_files_prefix']}{{contrast}}.txt")
    contrast_group_files_list = [os.path.join(config["contrast_group_files_prefix"], f"{contrast}.txt") for contrast in contrasts]
    contrasts_df_list = []
    for contrast, file in zip(contrasts, contrast_group_files_list):
        df = pd.read_csv(file, sep="\t", header=None, names=["sample", "Group"])
        df["ContrastName"] = contrast
        contrasts_df_list.append(df)
    contrasts_df = pd.concat(contrasts_df_list, ignore_index=True)
    contrasts_df = contrasts_df.merge(samples[["sample", "STARGenomeName"]], on="sample", how="left")
else:
    contrasts = []
    contrasts_df = pd.DataFrame(columns=["sample", "Group", "ContrastName", "STARGenomeName"])

# Exclude 1vs1 contrasts for differential expression targets
try:
    contrasts_de = [
        c for c in contrasts
        if contrasts_df[contrasts_df["ContrastName"] == c]
              .groupby("Group")["sample"].nunique().max() > 1
    ]
except NameError:
    contrasts_de = contrasts


# Input functions and other functions for the snakemake

def min_samples_per_group_for_contrast(wildcards):
    df = contrasts_df[contrasts_df['ContrastName'] == wildcards.contrast]
    group_counts = df.groupby('Group')['sample'].nunique()
    return min(group_counts.min(), 2)

def FillGenomeNameInFormattedString(FormattedString):
    def InputFunctionToReturn(wildcards):
        return FormattedString.format(GenomeName = samples[samples['sample']==wildcards.sample]['STARGenomeName'].iloc[0])
    return InputFunctionToReturn

def UsePairedEndFeatureCountsIfMixingSingleAndPairedReads(wildcards):
    is_single = samples.loc[samples['STARGenomeName']==wildcards.GenomeName]['Library_Layout'] == 'SINGLE'
    if is_single.any() and not is_single.all():
        return "-p"
    else:
        return ""

def GetProvidedBams(wildcards):
    return samples.loc[(samples['sample']==wildcards.sample) & samples['bam'].notna()]['bam'].tolist()

def GetDownloadLinks(wildcards):
    """
    ENA serves the same host and path over https as over ftp, and compute nodes can reach 443 but
    not 21, so rewriting the scheme is what lets rule DownloadFromAccession run off the login node
    """
    links = samples.loc[samples['sample']==wildcards.sample][wildcards.Read + '_link'].tolist()
    return [link.replace("ftp://", "https://", 1) if link.startswith("ftp://") else link for link in links]

def ExpandAllSamplesInFormatStringFromGenomeNameWildcard(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[samples['STARGenomeName']==wildcards.GenomeName]['sample'].unique())
    return InputFunctionToReturn

def ExpandAllSTARSamplesInFormatStringFromGenomeNameWildcard(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples_ForSTAR_df.loc[(samples_ForSTAR_df['STARGenomeName']==wildcards.GenomeName)]['sample'].unique())
    return InputFunctionToReturn

def ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[(samples['STARGenomeName']==wildcards.GenomeName) & (samples['Strandedness']==wildcards.Strandedness)]['sample'].unique())
    return InputFunctionToReturn

def ExpandAllSTARSamplesInFormatStringFromGenomeNameAndStrandWildcards(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples_ForSTAR_df.loc[(samples_ForSTAR_df['STARGenomeName']==wildcards.GenomeName) & (samples_ForSTAR_df['Strandedness']==wildcards.Strandedness)]['sample'].unique())
    return InputFunctionToReturn

def GetIndexSuffix(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] == True:
        return 'csi'
    else:
        return 'tbi'

def GetIndexingParamsFromGenomeName(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] == True:
        return '--csi'
    else:
        return ''

def GetIndexingParamsFromSampleName(wildcards):
    GenomeName = samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0]
    if STAR_genomes.loc[GenomeName]['ChromLargerThan512Mbp'] == True:
        return '--csi'
    else:
        return ''
    
def get_filled_path_from_contrast(wildcards, formatted_string):
    # Find the first row where 'ContrastName' matches wildcards.contrast
    row = contrasts_df[contrasts_df['ContrastName'] == wildcards.contrast].iloc[0]
    # Fill the formatted string with the 'STARGenomeName' from the found row
    return formatted_string.format(GenomeName=row['STARGenomeName'])

def GetMemForSuccessiveAttempts(*args, max_mb=48000):
    def ReturnMemMb(wildcards, attempt):
        i = int(attempt) - 1
        try:
            mem_mb = args[i]
        except IndexError:
            mem_mb = max_mb
        return mem_mb
    return ReturnMemMb

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 62000
