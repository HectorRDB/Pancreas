#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "segerstolpe"
loc="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_filt.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe_SC3.csv"
Rscript --verbose  3-sc3.R -n 32 -l $loc -o $out > 3-segerstolpe.out 2>&1
