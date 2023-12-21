import pandas as pd
import os


try:
    samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
    STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0)
except (NameError, KeyError) as NameOrKeyError:
    samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
    STAR_genomes = pd.read_csv("module_workflows/snakemake-workflow_rna-seq/config/STAR_Genome_List.tsv",sep='\t', index_col=0)


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

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 62000
