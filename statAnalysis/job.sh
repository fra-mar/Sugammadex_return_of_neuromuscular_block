#!/bin/bash

#SBATCH -A uppmaxPROJECTNUMBER
#SBATCH -t 06:00:00
#SBATCH -p pelle
#SBATCH -c 8 
 
#module load R/4.5.1 #Uncomment if Uppmax environment
date
Rscript regrModel_251207.R

echo ended regrModel

date
