require(biomaRt)
add_refseq <- function(x, y, z) {
    ensembl <- useMart("ensembl")
    ensembl <- useDataset(z, mart=ensembl)
    attr <- c("ensembl_gene_id", "entrezgene_id")
    all.entrezgene <- getBM(attributes = attr, values = "*", mart = ensembl) %>%
        unique %>%
        filter(!duplicated(ensembl_gene_id))
    left_join(x, all.entrezgene, by = setNames("ensembl_gene_id", y))
}
fetch_refseq <- function(specie) {
    ensembl <- useMart("ensembl")
    ensembl <- useDataset(specie, mart=ensembl)
    attr <- c("ensembl_gene_id", "entrezgene_id")
    getBM(attributes = attr, values = "*", mart = ensembl) %>%
        unique %>%
        filter(!duplicated(ensembl_gene_id))
}
