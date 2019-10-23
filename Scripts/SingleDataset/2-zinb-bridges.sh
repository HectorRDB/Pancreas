#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p RM
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

timestamp=$(date +"%Y%m%d-%H%M%S")
basename=2-zbinb_Pancreas-baron_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 10 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

echo "baron"
loc="/pylon5/ib5phhp/hectorrb/Pancreas/ProcessedData/baron_filt.rds"
meta="/home/hectorrb/Pancreas/Data/Baron/baron_meta.csv"
out="/pylon5/ib5phhp/Pancreas/ProcessedData/baron_zinbW.rds"
plot="/home/hectorrb/Pancreas/Figures/EDA/baron"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  Rscript --verbose  2-zinb.R -n 10 -l $loc -o $out -p $plot -m $meta > 2-baron.out 2>&1
