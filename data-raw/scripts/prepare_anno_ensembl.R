library(tidyverse)
library(rtracklayer)
source("scripts/add_refseq.R")

argv <- commandArgs(trailingOnly = TRUE)
org <- argv[1]
input <- argv[2]
release <- argv[3]
var <- paste0(org, ".Ensembl", release)
output <- paste0(var, ".RData")

anno <- import(input) %>%
    as.data.frame %>%
    filter(type == "transcript")

stopifnot(org %in% c("Hs", "Mm", "Rn", "Bt"))
if (org == "Hs") ensembl <- "hsapiens_gene_ensembl"
if (org == "Mm") ensembl <- "mmusculus_gene_ensembl"
if (org == "Rn") ensembl <- "rnorvegicus_gene_ensembl"
if (org == "Bt") ensembl <- "btaurus_gene_ensembl"
tx_2_all <- add_refseq(anno, "gene_id", ensembl) %>%
    dplyr::select(id = transcript_id,
                  ensembl_gene = gene_id,
                  symbol = gene_name,
                  entrez_id = entrezgene_id,
                  transcript_type = transcript_biotype)

# Save results
assign(var, tx_2_all)
save_cmd <- paste0("save(", var, ', file = "', output, '", compress = "xz")')
eval(parse(text = save_cmd))
