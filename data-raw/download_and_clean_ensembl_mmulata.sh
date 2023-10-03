#!/bin/bash

# Note: does not work with Ensembl release < 97(?)
#ftp://ftp.ensembl.org/pub/release-101/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_FASTA=Macaca_mulatta.Mmul_10.cdna.all.fa.gz

wget ${URL_FTP}/fasta/macaca_mulatta/cdna/${FILE_FASTA}
Rscript scripts/prepare_anno_ensembl_fasta.R Mmu ${FILE_FASTA} ${RELEASE}

#rm -rf ${FILE_GTF}
