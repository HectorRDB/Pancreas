#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/baron_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/baron_comp3"
Rscript --verbose  Dune.R -n 32 -l $loc -o $out -S "1.5.50" -C "10" -m "k_25" -r $rsec > baron3.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/segerstolpe_comp3"
Rscript --verbose  Dune.R -n 32 -l $loc -o $out -S "1.5.50" -C "10" -m "k_25" -r $rsec > segerstolpe3.out 2>&1
