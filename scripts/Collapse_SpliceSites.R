#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Collapse_Juncsfiles
# @created     : Tuesday Nov 28, 2023 21:12:39 CST
#
# @description : Pre-processing for re-implementing (https://doi.org/10.1093/nargab/lqab041)... Given input of regtools annotate (with splice donors and acceptors added) as I have done in my custom script, for each donor over a threshold, output a tsv file. Do same for acceptors.
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/GRCh38_GencodeRelease44Comprehensive.uniq.annotated.with_ss_scores.tsv.gz scratch/test.donors.bed scratch/test.acceptors.bed 10", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

fn_in <- args[1]
fn_out_donors <- args[2]
fn_out_acceptors <- args[3]
Threshold <- as.numeric(args[4])

library(tidyverse)

dat <- read_tsv(fn_in) %>%
    mutate(DonorPos_OneBased = ifelse(strand == "+", start+1, end-1)) %>%
    mutate(AcceptorPos_OneBased = ifelse(strand == "+", end-1, start+1))

Donors <- dat %>%
    group_by(chrom, DonorPos_OneBased, strand) %>%
    mutate(score = sum(score)) %>%
    distinct(chrom, DonorPos_OneBased, strand, score, .keep_all=T) %>%
    ungroup() %>%
    mutate(name = str_glue("{chrom}:{DonorPos_OneBased}:{strand}:D")) %>%
    mutate(start = DonorPos_OneBased-1) %>%
    dplyr::select(`#Chrom`=chrom, start, stop=DonorPos_OneBased, name, TotalJuncCounts=score, strand, Is_annotated=known_donor, SS_seq=`5ss_seq`, SS_score=`5ss_score`)

Donors %>%
    filter(TotalJuncCounts >= Threshold) %>%
    arrange(`#Chrom`, start) %>%
    write_tsv(fn_out_donors, col_names=T)

Acceptors <- dat %>%
    group_by(chrom, AcceptorPos_OneBased, strand) %>%
    mutate(score = sum(score)) %>%
    distinct(chrom, AcceptorPos_OneBased, strand, score, .keep_all=T) %>%
    ungroup() %>%
    mutate(name = str_glue("{chrom}:{AcceptorPos_OneBased}:{strand}:A")) %>%
    mutate(start = AcceptorPos_OneBased-1) %>%
    dplyr::select(`#Chrom`=chrom, start, stop=AcceptorPos_OneBased, name, TotalJuncCounts=score, strand, Is_annotated=known_acceptor, SS_seq=`3ss_seq`, SS_score=`3ss_score`)

Acceptors %>%
    filter(TotalJuncCounts >= Threshold) %>%
    arrange(`#Chrom`, start) %>%
    write_tsv(fn_out_acceptors, col_names=T)
