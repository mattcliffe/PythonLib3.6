# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:56:41 2018

@author: matthew
np.array([0.46,0.66,0.915,1.195,1.465,1.83,2.535,3.5,4.5,5.75,7.25,9,11,13,15,16.75])
"""

import numpy as np
import os
from dateutil.parser import parse

def CalcRC(PM,D):
    y = np.abs(np.log(D/PM)/(np.sqrt(2)*np.log(1.5)))
    G = 0.5*(1+0.14112821*y+0.08864027*y**2+0.02743349*y**3-0.00039446*y**4+0.00328975*y**5)**-8;
    Er = np.zeros([len(D)])
    
    for lmn in np.linspace(0,len(D)-1,len(D)):
        
        if D[lmn] <= PM:
            Er[lmn] = 100*(1-G[lmn])
    
        if D[lmn] >= PM:
            Er[lmn] = 100*G[lmn]

    return Er
            
            
def LoadFile(fname,path):
    os.chdir(path)
    
    
        
    #Load constants from file
    dstore = np.genfromtxt(fname,delimiter=',',skip_header=12,max_rows=4)
    volume = dstore[1,1:17]
    density = dstore[2,1:17]
    
    #load datetime 
    del dstore
    dstore = np.genfromtxt(fname,dtype=str,delimiter=',',skip_header=18,skip_footer=1,usecols=0)
    dt=[]
    for lmn in np.linspace(0,len(dstore)-1,len(dstore)):
        lmn = int(lmn)
        dt.append(parse(dstore[lmn]).time())
    #load data
    del dstore
    bindata = np.genfromtxt(fname,delimiter=',',skip_header=18,skip_footer=1,usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))
    sfr = np.genfromtxt(fname,delimiter=',',skip_header=18,skip_footer=1,usecols=(22))
    sampletime = np.genfromtxt(fname,delimiter=',',skip_header=18,skip_footer=1,usecols=(21))
    #Load precalculated PMr
    pmr = np.genfromtxt(fname,delimiter=',',skip_header=18,skip_footer=1,usecols=(24))
    
    
    
    return bindata, sfr, sampletime, dt, volume, density, fname, pmr
    
def CalcPM(bindata,volume,density,sfr,sampletime,RC):
    totalmass = np.sum(bindata*volume*density*RC/100,1)
    pm = totalmass/(sfr*sampletime)
    
    return pm


def running_mean(x, N):
    if len(x) >=N:
        cumsum = np.cumsum(np.insert(x, 0, 0)) 
        avg = np.zeros(len(x))
        avg[N-2:len(x)-1]=(cumsum[N:] - cumsum[:-N]) / float(N)
        avg[0:N] = np.cumsum(x[0:N])/N

    else:
        avg = np.cumsum(x)/N
    return avg

