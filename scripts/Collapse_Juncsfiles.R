#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Collapse_Juncsfiles
# @created     : Tuesday Nov 28, 2023 21:12:39 CST
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
fns_in <- args[-1]

library(tidyverse)
library(data.table)

dat <- fns_in %>%
    lapply(fread) %>%
    bind_rows()

dat %>%
    separate(V11, into=c("bSize1", "bSize2"), sep=",", convert=T, remove=F) %>%
    mutate(JuncStart = V2 + bSize1) %>%
    mutate(JuncEnd = V3 - bSize2) %>%
    distinct(JuncStart, JuncEnd, .keep_all=T) %>%
    select(-c("bSize1", "bSize2", "JuncStart", "JuncEnd")) %>%
    arrange(V1, V2, V3) %>%
    write_tsv(fn_out, col_names=F)
