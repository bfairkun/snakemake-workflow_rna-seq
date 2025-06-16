#!/usr/bin/env Rscript

######################################################################
# @author      : Your Name (your@email.com)
# @file        : SpliceER_FilterAndCalculateSSE.R
# @created     : 2025-06-09 11:11
#
# @description : 
######################################################################

# Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/test.SSE.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Alpha.sorted.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta1.sorted.bed.gz rna-seq/SplicingAnalysis/SplisER_Quantifications/GRCh38_GencodeRelease44Comprehensive/Donors.Beta2.sorted.bed.gz ", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_out <- args[1]
f_donor_alpha <- args[2]
f_donor_beta1 <- args[3]
f_donor_beta2 <- args[4]

library(tidyverse)
library(data.table)

bed6 <- fread(f_donor_alpha, select=c(1:6))
alpha <- fread(f_donor_alpha, drop=c(1,2,3,5,6)) %>%
    column_to_rownames("name") %>% as.matrix()
beta1 <- fread(f_donor_beta1, drop=c(1,2,3,5,6)) %>%
    column_to_rownames("name") %>% as.matrix()
beta2 <- fread(f_donor_beta2, drop=c(1,2,3,5,6)) %>%
    column_to_rownames("name") %>% as.matrix()

sum.counts <- alpha + beta1 + beta2

## Manually inspected alpha, beta1, and beta2 entries along with a bam in IGV to confirm they are correctly counting what I intend
# alpha["chr20:4691990:+:D", "X11_LN229_polyA_CP3PlusBPN_3160_1"]
# beta1["chr20:4691990:+:D", "X11_LN229_polyA_CP3PlusBPN_3160_1"]
# beta2["chr20:4691990:+:D", "X11_LN229_polyA_CP3PlusBPN_3160_1"]
# sum.counts["chr20:4691990:+:D", "X11_LN229_polyA_CP3PlusBPN_3160_1"]

# sum.counts[3, c("X2_LCL_polyA_Risdiplam_3160_1","X2_LCL_polyA_Risdiplam_316_1")]

SSE = alpha / sum.counts

## Filter out entries if SSE is not above some threshold in at least one sample
# SSE.filtered = SSE[apply(SSE, 1, function(x) any(x > 0.5, na.rm=TRUE)), , drop=FALSE]
# dim(SSE.filtered)
## But actually I realize this only filters out about half of entries... Mine as well just keep it all


## Write out SSE matrix as bed6+, with 4 sigfigs
SSE_out <- cbind(bed6, round(SSE, 4))
write_tsv(SSE_out, file = f_out)
