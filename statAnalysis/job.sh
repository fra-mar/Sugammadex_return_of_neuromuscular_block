#!/bin/bash

#SBATCH -A uppmax2025-2-378
#SBATCH -t 06:00:00
#SBATCH -p pelle
#SBATCH -c 8 
 
module load R/4.5.1
date
Rscript regrModel_251207.R

echo ended regrModel

date
#Rscript sensitivityAnalysis_RNBvsPKPDROCSUG.R
  
date
echo end
