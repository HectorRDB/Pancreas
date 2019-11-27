#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Garbage/Garb3_baron"
Rscript --verbose  garbage.R -n 20 -l $loc -o $out -S "1.5.50" -C "10" -m "k_25" > baron.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Garbage/Garb3_segerstolpe"
Rscript --verbose  garbage.R -n 20 -l $loc -o $out -S "1.5.50" -C "10" -m "k_25" > segerstolpe.out 2>&1
