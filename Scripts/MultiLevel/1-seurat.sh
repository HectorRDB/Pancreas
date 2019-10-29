#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_filt.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/MultiLevel/baron_Seurat.csv"
Rscript --verbose  1-seurat.R -l $loc -o $out > 1-baron.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_filt.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/MultiLevel/segerstolpe_Seurat.csv"
Rscript --verbose  1-seurat.R -l $loc -o $out > 1-segerstolpe.out 2>&1
