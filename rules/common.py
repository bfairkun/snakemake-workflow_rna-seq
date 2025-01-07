import pandas as pd
import numpy as np
import os
import requests
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.32")

# general functions
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

def populate_Library_Layout(row):
    if not pd.isna(row['R1']) and not pd.isna(row['R2']):
        return "PAIRED"
    elif not pd.isna(row['R1']) and pd.isna(row['R2']):
        return "SINGLE"
    elif row['R1_link'] and row['R2_link']:
        return "PAIRED"
    elif row['R1_link'] and not row['R2_link']:
        return "SINGLE"
    else:
        return ""

def wildcard_constraints_from_list(l):
    if len(l) == 0:
        return "DUMMY"
    else:
        return "|".join(l)

# read and process samples files
try:
    samples = pd.read_csv(config["samples"],sep='\t')
    STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0, converters={"ChromLargerThan512Mbp": lambda x: x.lower() in ["yes", "true", "y", "True", "T"]})
except (NameError, KeyError) as NameOrKeyError:
    samples = pd.read_csv("config/samples.tsv",sep='\t')
    STAR_genomes = pd.read_csv("config/STAR_Genome_List.tsv",sep='\t', index_col=0, converters={"ChromLargerThan512Mbp": lambda x: x.lower() in ["yes", "true", "y", "True", "T"]})

STAR_genomes['GenomeName'] = STAR_genomes.index

# Add column for R1_link and R2_link if not existing, and fill links if SRA_accession provided
if 'R1_link' not in samples.columns:
    samples['R1_link'] = pd.NA
if 'R2_link' not in samples.columns:
    samples['R2_link'] = pd.NA
samples[["R1_link", "R2_link"]] = samples.apply(
    lambda row: pd.Series(fetch_ena_links(row["SRA_accession"])) 
    if pd.isna(row["R1_link"]) or row["R1_link"] == "" 
    else pd.Series([row["R1_link"], row["R2_link"]]),
    axis=1
)

# Add column for SE vs PE if contains any NAs/blanks or column doesn't already exist
if 'R1' not in samples.columns:
    samples['R1'] = pd.NA
if 'R2' not in samples.columns:
    samples['R2'] = pd.NA
if 'Library_Layout' not in samples.columns:
    samples['Library_Layout'] = pd.NA
samples['Library_Layout'] = samples.apply(lambda row: populate_Library_Layout(row) if pd.isna(row['Library_Layout']) or row['Library_Layout'] == "" else row['Library_Layout'], axis=1)

# validate
validate(samples, "../schemas/samples.schema.yaml")
validate(STAR_genomes, "../schemas/STAR_Genome_List.schema.yaml")

AllSamples = samples['sample'].unique()

samples_SingleEnd = samples[samples['Library_Layout']=='SINGLE']['sample']
samples_PairedEnd = samples[samples['Library_Layout']=='PAIRED']['sample']

samples_from_links = samples.loc[samples["R1"].isna() & samples["R1_link"].notna()]['sample']
samples_from_local = samples.loc[samples["R1"].notna()]['sample']

samples_ForSTAR = samples.loc[samples['Aligner']=='STAR']['sample']

SingleEndSamples_wildcards_regex = wildcard_constraints_from_list(samples_SingleEnd)
PairedEndSamples_wildcards_regex = wildcard_constraints_from_list(samples_PairedEnd)


# Input functions and other functions for the snakemake

def FillGenomeNameInFormattedString(FormattedString):
    def InputFunctionToReturn(wildcards):
        return FormattedString.format(GenomeName = samples[samples['sample']==wildcards.sample]['STARGenomeName'].iloc[0])
    return InputFunctionToReturn

def UsePairedEndFeatureCountsIfMixingSingleAndPairedReads(wildcards):
    R2_isna = samples.loc[samples['STARGenomeName']==wildcards.GenomeName]['R2'].isna()
    if R2_isna.any() and not R2_isna.all():
        return "-p"
    else:
        return ""

def ExpandAllSamplesInFormatStringFromGenomeNameWildcard(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[samples['STARGenomeName']==wildcards.GenomeName]['sample'].unique())
    return InputFunctionToReturn

def ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[(samples['STARGenomeName']==wildcards.GenomeName) & (samples['Strandedness']==wildcards.Strandedness)]['sample'].unique())
    return InputFunctionToReturn

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
