#!/usr/bin/env Rscript

if(interactive()){
  args <- scan(text=
                 "rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/clustering/leafcutter_perind_numers.counts.gz rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/juncTableBeds/PSI_ByMax.bed rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/juncTableBeds/PSI.bed rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/juncTableBeds/JuncCounts.bed rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/juncTableBeds/PSIDenom.bed rna-seq/SplicingAnalysis/leafcutter/SalpingoecaRosetta_ensemblv_59/juncTableBeds/PSIByMaxDenom.bed", what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}


leafcutter_count_numers <- args[1]
PSIOut <- args[2]
PSIOutByMax <- args[3]
JuncOut <- args[4]
PSIDenomOut <- args[5]
PSIByMaxDenom <- args[6]

library(tidyverse)

## Move this to different script
# Make PSI tables
Count.Table.mat <- read.table(leafcutter_count_numers, sep = ' ', nrows = Inf) %>%
  as.matrix()

ClusterMax.mat <- Count.Table.mat %>%
  as.data.frame() %>%
  rownames_to_column("junc") %>%
  mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
  group_by(cluster) %>%
  mutate(across(where(is.numeric), max)) %>%
  ungroup() %>%
  select(junc, everything(), -cluster) %>%
  column_to_rownames("junc") %>%
  as.matrix()


ClusterSum.mat <- Count.Table.mat %>%
  as.data.frame() %>%
  rownames_to_column("junc") %>%
  mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
  group_by(cluster) %>%
  mutate(across(where(is.numeric), sum)) %>%
  ungroup() %>%
  select(junc, everything(), -cluster) %>%
  column_to_rownames("junc") %>%
  as.matrix()

PSI.df <- (Count.Table.mat / as.numeric(ClusterSum.mat) * 100) %>%
  signif() %>%
  as.data.frame()

PSIByWithinClusterMax.df <- (Count.Table.mat / as.numeric(ClusterMax.mat) * 100) %>%
  signif() %>%
  as.data.frame()


PSI.df %>%
    rownames_to_column("junc") %>%
    separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
    mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
    arrange(`#Chrom`, start, end) %>%
    write_tsv(PSIOut)

PSIByWithinClusterMax.df %>%
    rownames_to_column("junc") %>%
    separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
    mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
    arrange(`#Chrom`, start, end) %>%
    write_tsv(PSIOutByMax)

Count.Table.mat %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
    mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
    arrange(`#Chrom`, start, end) %>%
    write_tsv(JuncOut)

ClusterSum.mat %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
    mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
    arrange(`#Chrom`, start, end) %>%
    write_tsv(PSIDenomOut)

ClusterMax.mat %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
    mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
    arrange(`#Chrom`, start, end) %>%
    write_tsv(PSIByMaxDenom)
