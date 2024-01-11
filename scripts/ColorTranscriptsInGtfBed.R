#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79
# @file        :
# @created     : 2024-01-11 09:37:28
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text="../../ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Reference.gtf.bed.gz ../code/scratch/test.bed ../code/scratch/ColorKey.png", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

bed_in <- args[1]
bed_out <- args[2]
ColorKeyPdf_out <- args[3]

bed <- read_tsv(bed_in, col_names = c("chrom", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts", "gene_id","transcript_id","gene_type","gene_name","transcript_type","transcript_support_level","tag","transcript_name"  ),
                col_types=cols(.default = "c", start = "i", stop="i"))

bed %>%
  count(tag) %>%
  print(n=Inf)

bed %>%
  count(transcript_type, name="transcript_count") %>%
  print(n=Inf)

bed %>%
  count(gene_type, name="gene_count") %>%
  print(n=Inf)

full_join(
  bed %>%
    count(transcript_type, name="transcript_count"),
  bed %>%
    count(gene_type, name="gene_count"),
  by=c("transcript_type"="gene_type")) %>%
  print(n=Inf)

bed_out.df <- bed %>%
  mutate(HexColor = case_when(
    transcript_type == "protein_coding" ~ "#377eb8", #blue
    transcript_type == "nonsense_mediated_decay" ~ "#e41a1c", #red
    transcript_type == "retained_intron" ~ "#984ea3", #purple
    transcript_type == "processed_transcript" ~ "#f781bf", #pink
    transcript_type == "protein_coding_CDS_not_defined" ~ "#ff7f00", #orange
    transcript_type == "non_stop_decay" ~ "#a65628", #brown
    transcript_type == "protein_coding_LoF" ~ "#ffff33", #yellow
    transcript_type == "lncRNA" ~ "#4daf4a", #green
    TRUE ~ "#999999", #gray
  )) %>%
  mutate(itemRgb = apply(col2rgb(HexColor), 2, paste, collapse=','))

bed_out.df %>%
  distinct(HexColor, transcript_type) %>%
  ggplot(aes(x=1, y=transcript_type)) +
  geom_raster(aes(fill=HexColor)) +
  scale_fill_identity() +
  scale_x_continuous(expand=c(0,0))
ggsave(ColorKeyPdf_out, height=7, width=4)

bed_out.df %>%
  unite(Extras, gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_support_level, tag, transcript_name, sep=";") %>%
  dplyr::select(-HexColor) %>%
  write_tsv(bed_out, col_names = F)

