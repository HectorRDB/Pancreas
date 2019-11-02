#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Garbage/baron_"
Rscript --verbose  garbage.R -n 32 -l $loc -o $out -S "1.2.50" -C "0" -m "k_45" > baron1.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Garbage/segerstolpe_"
Rscript --verbose  garbage.R -n 32 -l $loc -o $out -S "1.2.50" -C "0" -m "k_45" > segerstolpe1.out 2>&1
