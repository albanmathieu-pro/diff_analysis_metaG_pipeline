#!/bin/bash

set -e
# Note: works with Ensembl release < 82

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_GTF=Rattus_norvegicus.Rnor_5.0.${RELEASE}.gtf.gz

wget ${URL_FTP}/gtf/rattus_norvegicus/${FILE_GTF}
Rscript scripts/prepare_anno_ensembl_old.R Rn ${FILE_GTF} ${RELEASE}

rm -rf ${FILE_GTF}
