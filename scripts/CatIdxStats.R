#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CatIdxStats
# @created     : Friday Dec 08, 2023 16:24:44 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


fn_out <- args[1]
fn_in <- args[-1]

library(tidyverse)

fn_in %>%
    lapply(read_tsv, col_names=c("chrom", "len", "counts", "Col4")) %>%
    bind_rows(.id="sample") %>%
    write_tsv(fn_out)
