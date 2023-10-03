library(tidyverse)
library(Biostrings)
source("scripts/add_refseq.R")

argv <- commandArgs(trailingOnly = TRUE)
org <- argv[1]
input <- argv[2]
release <- argv[3]
var <- paste0(org, ".Ensembl", release)
output <- paste0(var, ".RData")

anno_raw <- readDNAStringSet(input) %>% names
id <- str_extract(anno_raw, "^ENS[^\\.]*")
ensembl_gene <- str_extract(anno_raw, "gene:[^\\.]*") %>% str_replace("gene:", "")
symbol <- str_extract(anno_raw, "gene_symbol:[^ ]*") %>% str_replace("gene_symbol:", "")
transcript_type <- str_extract(anno_raw, "transcript_biotype:.*gene_symbol:") %>%
    str_replace("transcript_biotype:", "") %>%
    str_replace(" gene_symbol:", "")
tx_2_all <- tibble(id = id,
                  ensembl_gene = ensembl_gene,
                  symbol = symbol,
                  transcript_type = transcript_type)

stopifnot(org %in% c("Hs", "Mm", "Rn", "Mmu"))
if (org == "Hs") ensembl <- "hsapiens_gene_ensembl"
if (org == "Mm") ensembl <- "mmusculus_gene_ensembl"
if (org == "Rn") ensembl <- "rnorvegicus_gene_ensembl"
if (org == "Mmu") ensembl <- "mmulatta_gene_ensembl"
entrez_id <- fetch_refseq(ensembl)
tx_2_all <- left_join(tx_2_all, entrez_id, by = c("ensembl_gene" = "ensembl_gene_id")) %>%
    dplyr::select(id:symbol, entrez_id = entrezgene_id, transcript_type)

# Save results
assign(var, tx_2_all)
save_cmd <- paste0("save(", var, ', file = "', output, '", compress = "xz")')
eval(parse(text = save_cmd))
