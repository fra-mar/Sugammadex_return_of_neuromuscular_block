# Sugammadex_return_of_neuromuscular_block

Yes, this is complex if you are not used to programming and/or Linux. Feel free to contact me at francisco.martinez.torrente@akademiska.se for details on installation and execution.<br>
The aim of the whole repository is to made the code reproducible for all researches interested in the subject.

This repository contains the files needed to perform PKPD simulations resembling the administration of rocuronium and then reversal with sugammadex in a variety of custom scenarios.<br>
The code was designed to be run in a *High Performance Computing* environment (UPPMAX) due to the massive amount of data involved at some points. This environment is based on Linux and requires a bash file (.sh) to start with. Therefore the sequence of execution is dictated by "job.sh" file.<br>

You will need this software installed in your computer (all free): <br>
Python3 (version 3.12.3 used in this code) with the libraries / modules:<br>
  -*scipy*<br>
  -*numpy*<br>
  -*matplotlib*<br>
  -*pandas*<br>
R (I used version 4.5.1) with the packages/libraries:<br>
  -*ggplot2*<br>
  -*patchwork*<br>
  -*brms*<br>
  -*tidybayes*<br>
  -*bayesplot*<br>
  -*dplyr*<br>
  
  
  
This are the steps to reproduce the code. Do it first. Then change the scenarios if you want to test different regimes.<br>
<br>
1. **Create a project folder** and download the code to it. **Don't** change the folder structure.
2. Edit *job.sh*. If you are using your laptop/desktop aim for a *nSub* (i.e. number of subjects or number of simulated cases) under 2000. This number can be handled by most 8Gb RAM. If you have more than that feel free to experiment with a higher number.
3. Declare *job.sh* as an executable file by *chmod -x job.sh* and run it by typing *./job.sh*
4. job.sh will create a *simulations* folder and call the python code   *simulatorMain_251108.py*
5. This will gather the scenarios contained in *scenarios.csv*. For each scenario will create :<br>
   -a subfolder in *simulations*<br>
   -will gather PKPD parameters from *paramsBuilder_251019.py*<br>
   -will run the PKPD model from *kleijnModel_251104.py*<br>
   -will save a bunch of generated simulation data in the respective scenario subfolder.<br>
<br>Then will run the scripts *tofAnalyzer.R* to assess for RNB events, among other calculations, *summarizer.R* for summaries, and *plotter.R* for figures depicting concentrations at the central, peripheral and effect site compartments besides TOF ratio. Will plot all RNB events and max 200 no RNB events. 200 was chosen to allow for visualization. Otherwise the plot is a mess.<br>

Then switch to the *statAnalysis* folder <br>
There you run *regrModel_251207.R* which run the regression models in *brms* and saves the outputs to .RDS files.<br>
Then you run *mainStats_260202.R* that output stats, plots after reading the .RDS generated in the previous rows.<br>
<br>
Yes, it's a complex structure and difficult to grasp if not used to linux or running Python or R from command line, but that's how I managed to run RAM demanding calculations.
