library(tidyverse)
library(rtracklayer)
source("scripts/add_refseq.R")

argv <- commandArgs(trailingOnly = TRUE)
input <- argv[1]
release <- argv[2]
var <- paste0("Hs.Ensembl", release)
output <- paste0(var, ".RData")

anno <- import(input) %>%
    as.data.frame %>%
    filter(type == "transcript")

tx_2_all <- add_refseq(anno, "gene_id") %>%
    dplyr::select(id = transcript_id,
                  ensembl_gene = gene_id,
                  symbol = gene_name,
                  entrez_id = entrezgene,
                  transcript_type = transcript_biotype)

# Save results
assign(var, tx_2_all)
save_cmd <- paste0("save(", var, ', file = "', output, '", compress = "xz")')
eval(parse(text = save_cmd))
