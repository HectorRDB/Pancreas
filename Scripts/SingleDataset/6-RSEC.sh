#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron_Rsec.csv"
sce="/scratch/users/singlecell/Pancreas/ProcessedData/baron_Rsec.rds"
Rscript --verbose 6-RSEC.R -l $loc -o $out -s $sce > 6-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_zinbW.rds"
# out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/muraro_Rsec.csv"
# Rscript --verbose  6-Rsec.R -l $loc -o $out -s $sce > 6-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_zinbW.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe_Rsec.csv"
sce="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"
Rscript --verbose 6-RSEC.R -l $loc -o $out -s $sce > 6-segerstolpe.out 2>&1
