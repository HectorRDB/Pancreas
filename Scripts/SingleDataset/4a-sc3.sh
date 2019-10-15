#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/baron_filt_norm.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron_SC3.csv"
Rscript --verbose  4-sc3.R -n 32 -l $loc -o $out > 4-baron.out 2>&1

# echo "muraro"
# loc="/scratch/users/singlecell/Pancreas/ProcessedData/muraro_filt.rds"
# out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/muraro_SC3.csv"
# Rscript --verbose  3-sc3.R -n 32 -l $loc -o $out > 3-muraro.out 2>&1
