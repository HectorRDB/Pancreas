#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/baron.csv"
  Rscript --verbose  5-monocle.R -l $loc -o $out -K 1 > 5-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_zinbW.rds"
# out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/muraro.csv"
# Rscript --verbose  5-monocle.R -l $loc -o $out -k 3 > 5-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/segerstolpecsv"
  Rscript --verbose  5-monocle.R -l $loc -o $out -K 2 > 5-segerstolpe.out 2>&1
