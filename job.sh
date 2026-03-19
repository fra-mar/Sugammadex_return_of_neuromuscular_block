#!/bin/bash

#SBATCH -A xxxxxxxxxxxx # Your Uppmax project
#SBATCH -t 05:00:00
#SBATCH -p pelle 
#SBATCH -c 16 
 
#module load Python/3.12.3-GCCcore-13.3.0 # Uncomment for Uppmax environment
#module load R/4.5.1                      # Uncomment for Uppmax environment

rm -rf simulations
rm -rf outputs
rm -rf plots

mkdir simulations

timeStart=$(date +%Y-%m-%d\ %H:%M:%S)

echo -------------------
date
echo Run simulatorMain_251108
n=$(cat scenarios.csv | wc -l)

for i in $(seq 0 $(($n - 2)))
do
  echo ------------------------------
	echo $i
	python simulatorMain_251108.py $i 10000 # Number of simulations. If not high performace enviroment, <1500 recommended.
	date
done

echo -------------------
date 
echo Run tofAnalyzer.R
Rscript tofAnalyzer.R

echo -------------------
date 
echo Run summarizer
Rscript summarizer.R

echo ------------------
date
echo Run plotter
Rscript plotter.R
mkdir plots
cp ./simulations/scenario_*/*.png ./plots
echo ------------------
echo START $timeStart
date
echo END 
