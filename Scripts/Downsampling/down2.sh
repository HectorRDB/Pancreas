#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Downsampling/baron_comp2"
Rscript --verbose  downsampling.R -n 20 -l $loc -o $out -S "2.1.30" -C "0" -m "k_20" > baron2.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Downsampling/segerstolpe_comp2"
Rscript --verbose  downsampling.R -n 20 -l $loc -o $out -S "2.1.30" -C "0" -m "k_20" > segerstolpe2.out 2>&1
