#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/baron.csv"
  Rscript --verbose  01-zinb_tsne.R -l $loc -o $out -K 10 > 1-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_zinbW.rds"
# out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/muraro.csv"
# Rscript --verbose  01-zinb_tsne.R -l $loc -o $out -k 30 > 1-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/tSNE/segerstolpecsv"
  Rscript --verbose  01-zinb_tsne.R -l $loc -o $out -K 20 > 1-segerstolpe.out 2>&1
