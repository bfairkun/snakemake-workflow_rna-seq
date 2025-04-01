#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : BasicDifferentialExpression_edgeR
# @created     : Tuesday Mar 11, 2025 09:52:48 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/GRCh38_GencodeRelease44Comprehensive/AllSamplesUnstrandedCounting.Counts.txt config/contrast_group_files/BPN.vs.BPN_plus_SM.txt scratch/DE.tsv.gz scratch/plots", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

featCountsIn <- args[1]
contrastFile <- args[2]
outFile <- args[3]
plotsDir <- ifelse(length(args) > 3, args[4], NULL)

library(edgeR)
library(tidyverse)

# Create plots directory if it doesn't exist
if (!is.null(plotsDir)) {
    dir.create(plotsDir, showWarnings = FALSE, recursive = TRUE)
}

#Read in the contrast file
contrast <- read_tsv(contrastFile, col_names=F) %>%
    dplyr::rename("sample"=1, "Group"=2)

#Read in the featureCounts file
featCounts <- read_tsv(featCountsIn, comment="#") %>%
    dplyr::rename_with(~ str_replace(., ".+/(.+?)/Aligned.sortedByCoord.out.bam", "\\1"), -c(1:6)) %>%
    dplyr::select(1:6, unique(contrast$sample))

# Create a named factor for groups
group_factor <- factor(contrast$Group, levels = unique(contrast$Group))
names(group_factor) <- contrast$sample


#Create the DGEList object
y <- featCounts %>%
    dplyr::select(-c(2:6)) %>%
    column_to_rownames("Geneid") %>%
    DGEList(group=group_factor[colnames(.)])

# Filter low-expression genes
# keep <- filterByExpr(y)
# y <- y[keep, , keep.lib.sizes=FALSE]

mean.cpm <- calcNormFactors(y) %>%
    cpm(log=T) %>%
    apply(1, mean)

y <- y[mean.cpm > 1, , keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

cpm <- cpm(y, log=TRUE, prior.count=0.1)

# Estimate dispersion
y <- estimateCommonDisp(y)
y <- estimateTrendedDisp(y)
# default tagwise disp estimate is to shrink towards the trend
y <- estimateTagwiseDisp(y)

# Plot Biological Coefficient of Variation (BCV)
if (!is.null(plotsDir)) {
    png(file.path(plotsDir, "BCV_plot.png"))
    plotBCV(y)
    dev.off()}

# Create the design matrix with an intercept term
design <- model.matrix(~group_factor)

# Fit the model
fit <- glmQLFit(y, design)

# Perform differential expression analysis
qlf <- glmQLFTest(fit, coef=2)

# Save results to a file
tests.qlf <- topTags(qlf, n=Inf) %>%
    as.data.frame() %>%
    rownames_to_column("GeneID") %>%
    as_tibble()

write_tsv(tests.qlf, outFile)

tests.qlf.toplot <- tests.qlf %>%
    mutate(logFC = case_when(
        logFC > 3 ~ 3,
        logFC < -3 ~ -3,
        TRUE ~ logFC
    )) %>%
    mutate(
        shape = case_when(
            logFC > 3 ~ 24,
            logFC < -3 ~ 25,
            TRUE ~ 21
        ) 
    ) %>%
    mutate(shape = as.factor(shape))

N.tests <- tests.qlf.toplot %>%
    mutate(Sig = FDR<0.05) %>%
    count(Sig) %>%
    mutate(label=if_else(Sig, str_glue("{n} genes FDR<0.05"), str_glue("{n} genes FDR>=0.05"))) %>%
    mutate(vjust = 2*row_number())

volcano.plot <- tests.qlf.toplot %>%
    ggplot(aes(x=logFC, y=-log10(FDR), color=FDR<0.05)) +
    geom_point(aes(shape=shape)) +
    theme_minimal() +
    ggtitle("Volcano Plot") +
    xlab("log2 Fold Change") +
    ylab("-log10 Adjusted.P") +
    xlim(-3.1, 3.1) +
    scale_shape_manual(values=c("21"=21, "24"=24, "25"=25)) +
    geom_text(data = N.tests, x=-Inf, y=Inf, aes(label=label, vjust=vjust, color=Sig), hjust=0) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
    scale_color_manual(values=c("black", "red"), labels=c("FDR >= 0.05", "FDR < 0.05"), name=NULL) +
    theme(legend.position="none")

if (!is.null(plotsDir)) {
    ggsave(file.path(plotsDir, "Volcano_plot.png"), plot = volcano.plot)
}

#MA plot
MA.plot <- tests.qlf.toplot %>%
    ggplot(aes(x=logCPM, y=logFC, color=FDR<0.05)) +
    geom_point(aes(shape=shape)) +
    geom_text(data = N.tests, x=-Inf, y=Inf, aes(label=label, vjust=vjust, color=Sig), hjust=0) +
    theme_minimal() +
    ggtitle("MA Plot") +
    xlab("Log CPM") +
    ylab("log2 Fold Change") +
    scale_color_manual(values=c("black", "red"), labels=c("FDR >= 0.05", "FDR < 0.05")) +
    ylim(c(-3.1,3.1)) +
    scale_shape_manual(values=c("21"=21, "24"=24, "25"=25)) +
    theme(legend.position="none")

if (!is.null(plotsDir)) {
    ggsave(file.path(plotsDir, "MA_plot.png"), plot = MA.plot)
}
