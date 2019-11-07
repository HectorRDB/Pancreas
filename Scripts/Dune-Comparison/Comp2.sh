#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/baron_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/baron_comp2"
Rscript --verbose  Dune.R -n 32 -l $loc -o $out -S "2.1.30" -C "0" -m "k_20" -r $rsec > baron2.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/segerstolpe_comp2"
Rscript --verbose  Dune.R -n 32 -l $loc -o $out -S "2.1.30" -C "0" -m "k_20" -r $rsec > segerstolpe2.out 2>&1
