library(tidyverse)
library(rtracklayer)
gtf <- import("ERCC92.gtf")
ERCC92 <- tibble(id = gtf$gene_id,
                 ensembl_gene = gtf$gene_id,
                 symbol = as.character(seqnames(gtf)),
                 entrez_id = NA,
                 transcript_type = gtf$type)
save(ERCC92, file = "ERCC92.RData", compress = "xz")
