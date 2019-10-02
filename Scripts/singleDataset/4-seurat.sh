#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_filt.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron_Seurat.csv"
Rscript --verbose  4-seurat.R -l $loc -o $out > 4-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_filt.rds"
# out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/muraro_Seurat.csv"
# Rscript --verbose  4-seurat.R -l $loc -o $out > 4-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_filt.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe_Seurat.csv"
Rscript --verbose  4-seurat.R -l $loc -o $out > 4-segerstolpe.out 2>&1
