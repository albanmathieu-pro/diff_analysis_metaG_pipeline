#!/bin/bash

## Found on this page
## https://www.thermofisher.com/order/catalog/product/4456739#/4456739
wget -nc https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip

Rscript scripts/format_ERCC92.R

rm ERCC92.zip
rm ERCC92.fa
rm ERCC92.gtf
