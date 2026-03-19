#!/usr/bin/env python
# coding: utf-8
# Build params for Kleijn model implementation in Python

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import shutil
import json
import sys

roc2microM = 1.887505
sug2microM = 0.499451
# Fraction of parent RocBromide which is rocuronium
rocBromide2roc = 1 #529.8/609.7 

# for reproducibility
rng = np.random.default_rng(seed=1234)

# For same format as R
def rnorm(n, mean, sd):
    return rng.normal(mean,sd,n)


def buildParams(AGE, BW, CR, RAC, SEV, 
               bolusROC, i_ROC, durInf,
               bolusSUG, nSub):
    
      d0ROC  = bolusROC * BW            #induction dose ROC in mg (mgKg-1 *Kg)
      d0ROC  = d0ROC * rocBromide2roc * roc2microM #conversion to micromoles
      infROC = i_ROC * BW / 60                           # ROC mg per minute
      infROC = infROC * rocBromide2roc * roc2microM #conversion to micromoles
      d0SUG  = bolusSUG * BW            # Reversal dose SUG in mg (mgKg-1 *Kg)
      d0SUG  = d0SUG * sug2microM  # Conversion mg to micromoles
      tSUG   = durInf                # timestamp when SUG administered (hours)
      
      forEvents = pd.DataFrame(index=[0],
                               data={'d0ROC':d0ROC, 
                                     'infROC':infROC, 
                                     'durInf':durInf, 
                                     'd0SUG':d0SUG, 
                                     'tSUG':tSUG})
      
      CLage   = 1 + -0.00678 * (AGE - 43)
      tvCLro  = CLage * 0.269 *(BW/70)**0.75
      V1cr    = np.exp(-0.00143 * (CR-119))
      tvV1ro  = V1cr * 4.73 * (BW/70)**1
      Q2roRac = 1 if RAC=='nonAsian' else 1-0.212
      Q2ro    = Q2roRac * 0.279 * (BW/70)**0.75
      V2roAge = np.exp(0.00613 * (AGE-43))
      tvV2ro  = V2roAge * 6.76 * (BW/70)**1
      
      # for SUG
      CLsuBW  = 1 + 0.00378 * (BW-74.5)
      REN     = (2 * CR/(CR+119))**1.29
      tvCLsu  = CLsuBW * (REN * 0.093)
      V1suBW  = 1 + -0.00354 *(BW- 74.5)
      V1suRAC = 1 if RAC=='nonAsian' else 1-0.16
      tvV1su  = V1suBW * V1suRAC * 4.42 * (BW/70)**1
      tvQ2su  = 0.206 * (BW/70)**0.75
      V2suCR  = np.exp( -0.00305 * (CR-119))
      tvV2su  = V2suCR * 6.35 * (BW/70)**1
      
      # PD model ROC
      ke0Sev  = 1 if SEV==False else  1-0.567 # 1 no sevo,  1 +  -0.567 if sevo
      tvke0   = ke0Sev * 0.134 * (BW/70)**-0.25
      EC50sev = 1 if SEV==False else 1-0.395 # 1 no sevo,  1 +  -0.395 if sevo
      tvEC50  = EC50sev * 1.62 # microM
      #PD model SUG
      tvks    = np.exp(-3.43) * (BW/70)**-0.25
      tvHill  = 7.52
      tvE0    = 100 #1.04 * 100 #Arbitrary set to 100
      
      # eta
      eta_E0   = 0.01224571 # CV 11.1%
      eta_V1ro = 0.05509785 # CV% 23.8%
      eta_CLro = 0.09982362 # CV% 32.4%
      eta_V2ro = 0.09865367 # CV% 32.2%
      eta_ke0  = 0.16032217 #CV % 41.7%
      eta_EC50 = 0.06015486 # CV 24.9%
      eta_Hill = 0.15608110 # CV 41.1%
      eta_CLsu = 0.04895777 # CV 22.4%
      eta_ks   = 0.83273519 # CV 114%
      
      # Population estimates
      parsPop = pd.DataFrame(index=[0],
                             data = {
        'CLro' : tvCLro,
        'Q2ro' : Q2ro,
        'V1ro' : tvV1ro,
        'V2ro' : tvV2ro,
        'V1su' : tvV1su,
        'V2su' : tvV2su,
        'CLsu' : tvCLsu,
        'Q2su' : tvQ2su,
        'ke0' : tvke0,
        'EC50' : tvEC50,
        'Hill' : tvHill,
        'E0' : tvE0,
        'ks' : tvks}
      )
      # Log-normal random distributed parameters
      nSub = nSub - 1
      # To understand np.sqrt(eta_xxx) look at https://nmhelp.tingjieguo.com/V/4 4.1.1 
      parsRnd = pd.DataFrame(data = {
        'CLro' : tvCLro * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_CLro))),
        'Q2ro' : Q2ro,
        'V1ro' : tvV1ro * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_V1ro))),
        'V2ro' : tvV2ro * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_V2ro))),
        'V1su' : tvV1su,
        'V2su' : tvV2su,
        'CLsu' : tvCLsu * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_CLsu))),
        'Q2su' : tvQ2su,
        'ke0' : tvke0 * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_ke0))),
        'EC50' : tvEC50 * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_EC50))),
        'Hill' : tvHill * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_Hill))),
        'E0' : tvE0, # * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_E0))),
        'ks' : tvks * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_ks)))}
      )
      
      params = pd.concat((parsPop, parsRnd))
      params.reset_index(inplace=True, drop=True)
      return forEvents, params


if __name__ == "__main__":
    fe, params = buildParams(43, 74, 119,'nonAsian',True,
                             .6,.6,3.5,2,10000)
    print(fe)
    print(params.describe())

