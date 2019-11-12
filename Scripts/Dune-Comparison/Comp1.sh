#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

echo "baron"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/baron"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/baron_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/baron_comp1"
plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/ARI/baron_"
Rscript --verbose  Dune.R -n 16 -l $loc -o $out -S "1.2.50" -C "0" -m "k_50" -r $rsec -p $plot > baron1.out 2>&1

echo "segerstolpe"
loc="/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
rsec="/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/segerstolpe_comp1"
plot="/accounts/projects/epurdom/singlecell/Pancreas/Figures/ARI/segerstolpe_"
Rscript --verbose  Dune.R -n 16 -l $loc -o $out -S "1.2.50" -C "0" -m "k_50" -r $rsec -p $plot > segerstolpe1.out 2>&1
