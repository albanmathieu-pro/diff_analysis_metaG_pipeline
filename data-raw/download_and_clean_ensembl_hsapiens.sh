#!/bin/bash

# Note: does not work with Ensembl release < 97(?)

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_FASTA=Homo_sapiens.GRCh38.cdna.all.fa.gz

wget ${URL_FTP}/fasta/homo_sapiens/cdna/${FILE_FASTA}
#Rscript scripts/prepare_anno_ensembl_fasta.R Hs ${FILE_FASTA} ${RELEASE}
#
#rm -rf ${FILE_GTF}
