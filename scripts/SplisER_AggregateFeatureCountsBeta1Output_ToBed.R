#!/usr/bin/env Rscript

######################################################################
# @author      : Your Name (your@email.com)
# @file        : SplisER_AggregateFeatureCountsBeta1Output_ToBed.R
# @created     : 2025-06-06 09:24
#
# @description : 
######################################################################

# Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/test.out.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Alpha.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta1_U.Counts.txt rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta1_RF.Counts.txt rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta1_FR.Counts.txt", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_out <- args[1]
f_in_splicesites_by_samples_example_bed <- args[2] #Example input bed file in same format as output (same rows and columns, just different values for count table)
f_in <- args[-c(1, 2)]


library(data.table)
library(tidyverse)

bed6Plus_example <- fread(f_in_splicesites_by_samples_example_bed)

# Read all input featureCount files, keeping only columns 7:end for all except the first (which keeps 1:6)
dfs <- lapply(seq_along(f_in), function(i) {
  dat <- fread(f_in[i], header=TRUE, sep="\t", skip=1)
  dat <- dat %>%
    dplyr::rename_at(vars(-c(1:6)), ~str_replace(., "^.+/(.+?)/Aligned.sortedByCoord.out.bam", "\\1"))
  if (i == 1) {
    dat  # keep all columns for the first file
  } else {
    dat[, 7:ncol(dat), with=FALSE]  # keep only sample columns for others
  }
})

# Bind columns together
dat_merged <- bind_cols(dfs) %>%
  dplyr::select(name=Geneid, colnames(bed6Plus_example)[-c(1:6)])

bed6Plus_example %>%
  select(1:6) %>%
  inner_join(dat_merged, by="name") %>%
  write_tsv(f_out, col_names=TRUE)
