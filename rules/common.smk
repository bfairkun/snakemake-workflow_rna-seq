import pandas as pd
import numpy as np
import os


try:
    samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
    STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0)
except (NameError, KeyError) as NameOrKeyError:
    samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
    STAR_genomes = pd.read_csv("config/STAR_Genome_List.tsv",sep='\t', index_col=0)

samples_SingleEnd = samples[samples['R2'].isna()].index
samples_PairedEnd = samples[~samples['R2'].isna()].index

# # How to access values in samples.tsv

# print(samples)
# print( expand("Hello {sample}", sample=samples.index) )
# print( samples.at["A", "R1"] )

def FillGenomeNameInFormattedString(FormattedString):
    def InputFunctionToReturn(wildcards):
        return FormattedString.format(GenomeName = samples[samples.index==wildcards.sample]['STARGenomeName'].iloc[0])
    return InputFunctionToReturn

def ExpandAllSamplesInFormatStringFromGenomeNameWildcard(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[samples['STARGenomeName']==wildcards.GenomeName].index.unique())
    return InputFunctionToReturn

def ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=samples.loc[(samples['STARGenomeName']==wildcards.GenomeName) & (samples['Strandedness']==wildcards.Strandedness)].index.unique())
    return InputFunctionToReturn

def GetIndexingParamsFromGenomeName(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] == 'T':
        return '--csi'
    else:
        return ''

def GetIndexingParamsFromSampleName(wildcards):
    GenomeName = samples.loc[wildcards.sample]['STARGenomeName']
    if STAR_genomes.loc[GenomeName]['ChromLargerThan512Mbp'] == 'T':
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
