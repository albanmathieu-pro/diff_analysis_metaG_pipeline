#!/bin/bash

# Note: does not work with Ensembl release < 82

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_GTF=Bos_taurus.ARS-UCD1.2.${RELEASE}.chr.gtf.gz

wget ${URL_FTP}/gtf/bos_taurus/${FILE_GTF}
Rscript scripts/prepare_anno_ensembl.R Rn ${FILE_GTF} ${RELEASE}

rm -rf ${FILE_GTF}
