#!/usr/bin/env python
# coding: utf-8



from scipy.integrate import odeint
import numpy as np
import pandas as pd

def model(y, t, params, I):
    
    rocCentral, rocPerif, sugCentral, sugPerif, C1cx, C2cx, Ceff = y
    
    CLro, Q2ro, V1ro, V2ro = params[0:4]
    V1su, V2su, CLsu, Q2su = params[4:8]
    ke0, ks                = params[8], params[12]
    
    
    C1ro = rocCentral / V1ro
    #C2ro = rocPerif   / V2ro
    #SUG, PK      
    C1su = sugCentral / V1su
    #C2su = sugPerif / V2su
    #Complex, pK
    #C1cx = cpxCentral / V1su #complex same Kel and V1,2 as sug
    #C2cx = cpxPerif   / V2su
    K2 = 0.03404745 # in paper log(e) K2 = -3.38, ie exp(-3.38)
    K1 = 0.6090779 # in paper k1cx=k2cx/kd
    # Reparameterize Q2s
    k12ro = Q2ro / V1ro
    k21ro = Q2ro / V2ro

    k12su = Q2su / V1su
    k21su = Q2su / V2su
     
    ddt_rocCentral  = - (C1ro * C1su * K1) * V1su  + I  
    ddt_rocCentral +=   (C1cx * K2) * V1su - CLro * C1ro  
    ddt_rocCentral +=   k21ro * rocPerif - k12ro * rocCentral 
    ddt_rocPerif    =   k12ro * rocCentral - k21ro * rocPerif
  
    ddt_sugCentral  = - (C1ro * C1su * K1) * V1su  
    ddt_sugCentral +=   (C1cx * K2) * V1su - CLsu * C1su 
    ddt_sugCentral +=   k21su * sugPerif - k12su * sugCentral 
    ddt_sugPerif    =   k12su * sugCentral - k21su * sugPerif
  
    ddt_C1cx        =   (C1ro * C1su * K1) - (C1cx * K2) - (C1cx * CLsu/V1su)
    ddt_C1cx       +=   k21su * C2cx - k12su * C1cx 
    ddt_C2cx        =   k12su * C1cx - k21su * C2cx 

    Ceff_dt         =   ke0 * (C1ro - Ceff) - ks * C1su

    return (ddt_rocCentral, ddt_rocPerif,
            ddt_sugCentral, ddt_sugPerif,
            ddt_C1cx, ddt_C2cx, Ceff_dt)       

if __name__ == "__main__":
    from paramsBuilder_251019 import buildParams
    #import matplotlib.pyplot as plt
    
    fe, params = buildParams(43, 74, 119,'nonAsian',True,
                             .6,.6,3.5,2,100)
    print(fe)
    print(params.describe())

