#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Downsampling/baron_comp1"
Rscript --verbose  downsampling.R -n 20 -l $loc -o $out -S "1.2.50" -C "0" -m "k_50" > baron1.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Downsampling/segerstolpe_comp1"
Rscript --verbose  downsampling.R -n 20 -l $loc -o $out -S "1.2.50" -C "0" -m "k_50" > segerstolpe1.out 2>&1
