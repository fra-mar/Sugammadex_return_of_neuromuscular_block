#!/usr/bin/env python
# coding: utf-8

# Kleijn model for interaction ROC-SUG 
# Based on Britsh Journal Clinical Pharmacology 2011; 72: 415-433
# Results have been compared to same model implemented in R library RxODE2
# https://www.omnicalculator.com/chemistry/molarity
# SUG molar mass: 2002.2 g/mol https://pubchem.ncbi.nlm.nih.gov/compound/sugammadex
# ROC molar mass: 529.8 g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Rocuronium

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
from kleijnModel_251104 import model
from paramsBuilder_251019 import buildParams
import gc

args = sys.argv
if len(args) == 1:
    print('Script called without arguments for scenario number and nSubs')
else:
    # Which scenario (0...n scenarios)
    scenN = int(args[1])
    # N subjects
    nSub  = int(args[2])

# Load scenarios
df_scenarios = pd.read_csv('scenarios.csv')

# Follow up after SUG administration
follup = 60 * 3 # 60min * n hours. Units: minutes

# Select the scenario numbered when the script was called
s = df_scenarios.iloc[scenN,:]
# 'simulations' folder should exist / created previously by job.sh script
scenDir = os.path.join('simulations',f'scenario_{s.sceName}')
os.mkdir(scenDir)


# Load variables for events (fe), and pkpd parameters (params)
fe, params = buildParams(s.AGE, s.BW, s.CR,s.RAC,s.SEV,
                         s.d0ROC,s.infROC,s.infDur,
                         0,nSub)

timeline  = np.linspace(0,fe.loc[0,'durInf']*60,int(fe.loc[0,'durInf']*60)*6)
timeline_ = np.linspace(0, follup,follup*6) # follup defined in line 28-29
timeline_ = timeline_ + timeline[-1]
TIMELINE = np.hstack( (timeline, timeline_) ).reshape(-1,1)

dROC_tosave = s.d0ROC * s.BW + s.infROC * s.BW * s.infDur #to be used later

for i in range(nSub):
    #_----------------------------ROC phase ------------------------------------
    
    # Select parameters for the Subject[i]
    pars = params.iloc[i,:]
    
    # Load ROC and SUG doses
    d0ROC  = fe.loc[0,'d0ROC']  # Units mg (weight included)
    infROC = fe.loc[0,'infROC'] # Units mg (weight included) per minute
    
    # Initial values, ROC phase. AMOUNTS
    rocCentral, rocPerif = d0ROC, 0
    sugCentral, sugPerif = 0, 0
    C1cx, C2cx           = 0, 0
    Ceff                 = 0

    
    # Initial values, ROC phase. CONCENTRATIONS
    inits = (rocCentral, rocPerif,
             sugCentral, sugPerif, 
             C1cx, C2cx,
             Ceff)  
    
    
    # Run the diff equations
    r = odeint(model, inits, timeline, args=(np.array(pars),infROC))
    ROCcentral, ROCperif, _, _, _, _, C_ROCeff = r.T # Capital letters to human-read them as np.arrays
    # Get last values as initial values for SUG phase
    inits = r[-1,:]
    del r

    # Calculated TOFr
    Emax = pars.E0
    tofR = pars.E0 - ( (Emax * C_ROCeff**pars.Hill) / (pars.EC50**pars.Hill + C_ROCeff**pars.Hill)  )
    tofR = np.nan_to_num(x=tofR, nan=pars.E0) #when high TOF line above can RuntimeWarning
        
    # -------------------------SUG Phase---------------------------------------
    lastTOFr = tofR[-1]
    
    if "fixed" in s.d0SUG:
        d0SUG = float(s.d0SUG[5:]) / s.BW
    elif s.d0SUG == "var":
        d0SUG = 2.0 if lastTOFr > 10 else 4.0
    else:
        print('found nothing about SUG dose')
    
    dSUG_tosave = d0SUG # Spare for later comparison
    d0SUG = d0SUG * s.BW * .499451 #molar conversion



    # Get last values will be inits for SUG step
    rocCentral, rocPerif, sugCentral, sugPerif, C1cx, C2cx, Ceff = inits 
    sugCentral = d0SUG  # Correct the sugCentral with the actual dose SUG 
    inits = (rocCentral, rocPerif,
             sugCentral, sugPerif,
             C1cx, C2cx, 
             Ceff) 
    
    r = odeint(model, inits, timeline_, args=(np.array(pars),0))

    ROCcentral_s, ROCperif_s, _, _, _, _, C_ROCeff_s = r.T

    # Calculate TOFr for SUG phase
    tofR_ = pars.E0 - ( (Emax * C_ROCeff_s**pars.Hill) /
                       (pars.EC50**pars.Hill +C_ROCeff_s**pars.Hill)  )
    
    tofR_ = np.nan_to_num(x=tofR_, nan=pars.E0) #when high TOF line above can RuntimeWarning
        
      # Transform amounts to concentrations

    if i == 0:
        C1ro = np.hstack(( ROCcentral, ROCcentral_s)).reshape(-1,1) /pars.V1ro 
        C2ro = np.hstack(( ROCperif,   ROCperif_s)).reshape(-1,1)   /pars.V2ro 
        Ceff_ro = np.hstack(( C_ROCeff,   C_ROCeff_s)).reshape(-1,1)  # Risk confusion with Ceff!!!
        TOFr = np.hstack((tofR, tofR_)).reshape(-1,1) 

        forRegression = np.array([[int(i+1), dROC_tosave,
                                   dSUG_tosave,
                                   ROCcentral[-1],
                                   ROCcentral[-1]/pars.V1ro,
                                   ROCperif[-1],
                                   ROCperif[-1]/pars.V2ro]])
    
    else:
        C1ro = np.hstack(   (C1ro, 
                             np.hstack(( ROCcentral, ROCcentral_s)).reshape(-1,1) /pars.V1ro)   )
        C2ro = np.hstack(   (C2ro,
                             np.hstack(( ROCperif,   ROCperif_s)).reshape(-1,1)   /pars.V2ro)   )
        Ceff_ro = np.hstack((Ceff_ro,
                             np.hstack(( C_ROCeff,   C_ROCeff_s)).reshape(-1,1) ))  # Risk confusion with Ceff!!!
        TOFr = np.hstack((TOFr,
                          np.hstack((tofR, tofR_)).reshape(-1,1))  )
        
        forRegression = np.vstack( (forRegression,
                                   np.array([[int(i+1), 
                                               dROC_tosave, 
                                               dSUG_tosave,
                                               ROCcentral[-1],
                                               ROCcentral[-1]/pars.V1ro,
                                               ROCperif[-1],
                                               ROCperif[-1]/pars.V2ro]]) ) )
          

  
    #------------------End of for loop#

ids = (np.arange(0,nSub) +1).astype('str')
ids = ['id'+x for x in ids]
ids = ",".join( ids )
np.savetxt(os.path.join(scenDir,'C1ro.csv'), 
           C1ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'C2ro.csv'), 
           C2ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'Ceff_ro.csv'), 
           Ceff_ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'TOFr.csv'), 
           TOFr, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'timeline.csv'), 
           TIMELINE, delimiter=',',header='time', comments='')
np.savetxt(os.path.join(scenDir,'forRegression.csv'),
           forRegression, delimiter=',', header='id,dROC,dSUG,A1ro,C1ro,A2ro,C2ro', comments='')


# Save params
params_idx = np.array(params.index) +1
params.reset_index(drop=True, inplace=True)
params = pd.concat( (pd.DataFrame(data={'simId':params_idx}), params), axis=1)
params.to_csv(os.path.join(scenDir,'params.csv'), index=False)

print(f'\tDone with {scenDir}')
del C1ro, C2ro, Ceff_ro, TOFr, params, TIMELINE
gc.collect()


