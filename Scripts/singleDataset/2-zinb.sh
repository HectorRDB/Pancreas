#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_filt.rds"
out="/scratch/users/singlecell/Pancreas/ProcessedData/baron_zinbW.rds"
plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/EDA/baron"
Rscript --verbose  2-zinb.R -n 32 -l $loc -o $out -p plot > 2-baron.out 2>&1

echo "muraro"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_filt.rds"
out="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_zinbW.rds"
plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/EDA/muraro"
Rscript --verbose  2-zinb.R -n 32 -l $loc -o $out -p plot > 2-muraro.out 2>&1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_filt.rds"
out="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_zinbW.rds"
plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/EDA/segerstolpe"
Rscript --verbose  2-zinb.R -n 32 -l $loc -o $out -p plot > 2-segerstolpe.out 2>&1
# 
# echo "xin"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/xin_filt.rds"
# out="/scratch/users/singlecell/Pancreas/ProcessedData/xin_zinbW.rds"
# plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/EDA/xin"
# Rscript --verbose  2-zinb.R -n 32 -l $loc -o $out -p plot > 2-xin.out 2>&1
