#!/bin/bash
datasets="baron segerstolpe"

for dataset in $datasets
do
  echo $dataset
  Rscript --no-save --quiet render.R -d $dataset
done
