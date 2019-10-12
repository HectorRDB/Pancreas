#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/RawData/baron.rds"
out="/scratch/users/singlecell/Pancreas/ProcessedData/baron_filt.rds"
  meta="/accounts/projects/epurdom/singlecell/Pancreas/Data/Baron/baron_meta.csv"
Rscript --verbose  1-filtering.R -l $loc -o $out -c 5 -m $meta > 1-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/RawData/muraro.rds"
# out="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_filt.rds"
# Rscript --verbose  1-filtering.R -l $loc -o $out -c 5 > 1-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/RawData/segerstolpe.rds"
out="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_filt.rds"
meta="/accounts/projects/epurdom/singlecell/Pancreas/Data/Segerstolpe/segerstolpe_meta.csv"
Rscript --verbose  1-filtering.R -l $loc -o $out -c 5 -m $meta > 1-segerstolpe.out 2>&1
#
# echo "xin"
# loc="/scratch/users/singlecell/Pancreas/RawData/xin.rds"
# out="/scratch/users/singlecell/Pancreas/ProcessedData/xin_filt.rds"
# Rscript --verbose  1-filtering.R -l $loc -o $out -c 5 > 1-xin.out 2>&1
