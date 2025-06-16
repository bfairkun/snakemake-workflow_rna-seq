#!/usr/bin/env Rscript

######################################################################
# @author      : Your Name (your@email.com)
# @file        : SpliseER_DESeq2_SpliceSiteUsageContrast.R
# @created     : 2025-06-10 11:43
#
# @description : Use DESeq2 to compare Splice site usage (alpha vs beta1+beta2) between groups.
######################################################################

# Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "/project2/yangili1/bjf79/sm_splicingmodulators/code/config/Contrasts/X11_DMSO_v_CP3PlusBPN316.txt rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Alpha.sorted.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta1.sorted.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta2.sorted.bed.gz scratch/Donors.DESeq2.results.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(DESeq2)
library(tidyverse)

# --- Inputs ---
spliced_file <- "spliced.bed"
partial_unspliced_file <- "partial_unspliced.bed"
unspliced_file <- "unspliced.bed"
groups_file <- "groups.txt"

# --- Load group info ---
groups <- read_tsv(groups_file, col_names = c("sample", "group"))

# --- Helper function to read bed6+ files ---
read_bed6_plus <- function(file) {
  df <- read_tsv(file, col_names = FALSE)
  meta <- df[, 1:6]
  counts <- df[, -(1:6)]
  colnames(counts) <- make.unique(colnames(counts))
  rownames(counts) <- apply(meta, 1, paste, collapse = "_")
  counts
}

# --- Read all three files ---
spliced <- read_bed6_plus(spliced_file)
partial_unspliced <- read_bed6_plus(partial_unspliced_file)
unspliced <- read_bed6_plus(unspliced_file)

# --- Combine unspliced sources ---
total_unspliced <- partial_unspliced + unspliced

# --- Filter to samples of interest ---
samples <- groups$sample
spliced <- spliced[, samples, drop=FALSE]
total_unspliced <- total_unspliced[, samples, drop=FALSE]

# --- Stack spliced and unspliced into one matrix ---
all_counts <- rbind(spliced, total_unspliced)

# --- Prepare sample table ---
condition <- rep(groups$group, each = 2)
intron_status <- rep(c("spliced", "unspliced"), times = length(samples))
sample_id <- rep(samples, each = 2)

col_data <- data.frame(
  sample = factor(sample_id),
  group = factor(condition[seq(1, length(condition), by=2)]),
  intron_status = factor(intron_status)
)

# --- Build DESeq2 dataset ---
dds <- DESeqDataSetFromMatrix(
  countData = all_counts,
  colData = col_data,
  design = ~ sample + intron_status + group:intron_status
)

# --- Run DESeq2 ---
dds <- DESeq(dds)

# --- Results: interaction term tests if proportion differs by group ---
res <- results(dds, name = "groupB.intron_statusspliced")  # Assuming group A vs B

# --- Output ---
res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), file = "deseq2_intron_retention_results.csv")
