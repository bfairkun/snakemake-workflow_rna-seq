#!/usr/bin/env Rscript

######################################################################
# @author      : Your Name (your@email.com)
# @file        : leafcutter_tidyDifferentialSplicingResults.R
# @created     : 2025-06-09 12:05
#
# @description : 
######################################################################

# Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "rna_seq/SplicingAnalysis/ObservedJuncsAnnotations/GRCh38_GencodeRelease44Comprehensive.uniq.annotated.with_ss_scores.tsv.gz rna_seq/SplicingAnalysis/ClassifyJuncs/GRCh38_GencodeRelease44Comprehensive/Leaf2_junction_classifications.txt rna_seq/SplicingAnalysis/differential_splicing/Colombo.SM67_and_SMG7_vs_Scramble/leaf_effect_sizes.txt rna_seq/SplicingAnalysis/differential_splicing/Colombo.SM67_and_SMG7_vs_Scramble/leaf_cluster_significance.txt scratch/Results.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_junctions <- args[1]
f_classifications <- args[2]
f_effect_sizes <- args[3]
f_significance <- args[4]
f_out <- args[5]

library(tidyverse)
library(data.table)

junctions <- fread(f_junctions)
classifications <- fread(f_classifications)
effect_sizes <- fread(f_effect_sizes) %>%
    separate(intron, into=c("chrom", "start", "end", "cluster"), sep=":", remove=F, convert=T) %>%
    mutate(cluster = str_glue("{chrom}:{cluster}"))
significance <- fread(f_significance)

Joined <- effect_sizes %>%
    inner_join(significance, by=c("cluster")) %>%
    mutate(strand = str_replace(cluster, ".+?([+-])$", "\\1")) %>%
    inner_join(junctions, by=c("chrom", "start", "end", "strand")) %>%
    mutate(Intron_coord = str_glue("{chrom}:{start}-{end-1}")) %>%
    inner_join(
        classifications,
        by=c("Intron_coord")
    )

Joined %>%
    write_tsv(f_out, col_names=T)