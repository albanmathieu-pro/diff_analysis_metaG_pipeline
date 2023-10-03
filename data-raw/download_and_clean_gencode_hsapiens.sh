#!/bin/bash

RELEASE=$1
URL_FTP=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}
FILE_FA=gencode.v${RELEASE}.transcripts.fa.gz

#wget -nc ${URL_FTP}/${FILE_FA}

Rscript scripts/prepare_anno_gencode_hsapiens.R ${FILE_FA} ${RELEASE}
