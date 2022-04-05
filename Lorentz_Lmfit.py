#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code performs the lorentzian fit on the spectra
V:2022/04/01

"""

import numpy as np
from lmfit import minimize, Parameters

# OPTIONS
year='2015'
month='1503'
metodo= 'leastsq' 
# DATA
nventa=2**13
fm=256
df=float(fm)/float(nventa)
fajus_inf,fajus_sup= 6.35, 23.75
fcal_inf,fcal_sup= 6., 25.
nreso= 3
npara= nreso*3+2
def pofre(f,f0,df):
    return int(round((f-f0)/df))
fajus_inf_pos= pofre(fajus_inf,fcal_inf,df)
fajus_sup_pos= pofre(fajus_sup,fcal_inf,df) 
fre= np.arange(fcal_inf,fcal_sup+df,df)  
nf= len(fre)                                 
freajus= fre[fajus_inf_pos : fajus_sup_pos+1]  
nfajus= len(freajus)
        
# Paths
pathg="S_N_FD/"
pathmonth=pathg+year+"/"+month+"/"+"SR"+month+"_"
        
# Files (output from antropo.py) are loaded
rs0= np.genfromtxt(pathmonth+"mediaNA"+'_0').reshape(-1,nf)
rs0_ajus= rs0[:,fajus_inf_pos : fajus_sup_pos+1]
rs1=np.genfromtxt(pathmonth+"mediaNA"+'_1').reshape(-1,nf)
rs1_ajus= rs1[:,fajus_inf_pos : fajus_sup_pos+1]
nintervmes=len(rs0)
 
# Initial values
f1i, f2i, f3i= 8.011, 14.2, 20.63
s1i, s2i, s3i= 1.78, 1.94, 2.56
mi, ni= 0., 0.
para_ini0= np.array([0.,0.,0.,f1i,f2i,f3i,s1i,s2i,s3i,mi,ni])
para_ini1= np.array([0.,0.,0.,f1i,f2i,f3i,s1i,s2i,s3i,mi,ni])
f1i_pos= pofre(f1i,fajus_inf,df)
f2i_pos= pofre(f2i,fajus_inf,df)
f3i_pos= pofre(f3i,fajus_inf,df)

rs0_ajus_mediaMes= np.mean(rs0_ajus,0)
rs1_ajus_mediaMes= np.mean(rs1_ajus,0)


para_ini0[0]= rs0_ajus_mediaMes[f1i_pos] 
para_ini0[1]= rs0_ajus_mediaMes[f2i_pos] 
para_ini0[2]= rs0_ajus_mediaMes[f3i_pos]
para_ini1[0]= rs1_ajus_mediaMes[f1i_pos] 
para_ini1[1]= rs1_ajus_mediaMes[f2i_pos] 
para_ini1[2]= rs1_ajus_mediaMes[f3i_pos]
            
#params = Parameters()
params0 = Parameters()
params1 = Parameters()
    
def residual(params, x, data):
    v= params.valuesdict()     
    model = v['a1']/(1+((x-v['f1'])**2)/v['s1']**2)+\
                v['a2']/(1+((x-v['f2'])**2)/v['s2']**2)+\
                v['a3']/(1+((x-v['f3'])**2)/v['s3']**2)+v['m']*x+v['n']      
    return (data-model)
# A dictionary containing the parameters is created
def defpar(params,val):
    params.add('a1', value = val[0])
    params.add('a2', value = val[1])
    params.add('a3', value = val[2])    
    params.add('f1', value = val[3])
    params.add('f2', value = val[4])
    params.add('f3', value = val[5])    
    params.add('s1', value= val[6])
    params.add('s2', value= val[7])
    params.add('s3', value= val[8])    
    params.add('m', value= val[9])
    params.add('n', value= val[10])
# Output reading functions
def leepar(salida):
    val= np.zeros(11)
    val[0]=salida.params['a1'].value
    val[1]=salida.params['a2'].value
    val[2]=salida.params['a3'].value
    val[3]=salida.params['f1'].value
    val[4]=salida.params['f2'].value
    val[5]=salida.params['f3'].value    
    val[6]=salida.params['s1'].value
    val[7]=salida.params['s2'].value
    val[8]=salida.params['s3'].value
    val[9]=salida.params['m'].value
    val[10]=salida.params['n'].value
    return val
    
def leeerr(salida):
    val= np.zeros(11)
    val[0]=salida.params['a1'].stderr
    val[1]=salida.params['a2'].stderr
    val[2]=salida.params['a3'].stderr
    val[3]=salida.params['f1'].stderr
    val[4]=salida.params['f2'].stderr
    val[5]=salida.params['f3'].stderr   
    val[6]=salida.params['s1'].stderr
    val[7]=salida.params['s2'].stderr
    val[8]=salida.params['s3'].stderr
    val[9]=salida.params['m'].stderr
    val[10]=salida.params['n'].stderr
    return val
    
def leechi(salida):
    return salida.chisqr
    
        
# The parameter values for each 10-min interval are calculated
defpar(params0, para_ini0)
defpar(params1, para_ini1)

para_s0= np.zeros((nintervmes, npara),dtype=float)
para_s1= np.zeros((nintervmes, npara),dtype=float)
salerr0=np.zeros((nintervmes, npara),dtype=float)
salerr1=np.zeros((nintervmes, npara),dtype=float)
salchi0=np.zeros((nintervmes),dtype=float)
salchi1=np.zeros((nintervmes),dtype=float)
    
for i in np.arange(nintervmes):
   
# Sensor 0 
    out_s0 = minimize(residual, params0,\
            args=(freajus, rs0_ajus[i]), method = metodo)
    para_s0[i]= leepar(out_s0)
    salerr0[i]= leeerr(out_s0)
    salchi0[i]= leechi(out_s0)
# Sensor 1    
    out_s1 = minimize(residual, params1,\
            args=(freajus, rs1_ajus[i]), method = metodo)
    para_s1[i]= leepar(out_s1)
    salerr1[i]= leeerr(out_s1)
    salchi1[i]= leechi(out_s1)
        
salida0=np.zeros((nintervmes, npara+nreso),dtype=float)
salida1=np.zeros((nintervmes, npara+nreso),dtype=float)
salida0[:,nreso:]=para_s0[:,:]
salida1[:,nreso:]=para_s1[:,:]
# Global mode amplitudes
for i in np.arange(nintervmes):
    ampm0 = np.zeros(nreso)
    ampm1 = np.zeros(nreso)
    for k in np.arange(nreso):
        ampm0[k]=para_s0[i,nreso*3]*para_s0[i,nreso+k]+para_s0[i,nreso*3+1]
        ampm1[k]=para_s1[i,nreso*3]*para_s1[i,nreso+k]+para_s1[i,nreso*3+1]
        for l in np.arange(nreso):
            ampm0[k]+= para_s0[i,l]/(1+((para_s0[i,nreso+k]-\
                para_s0[i,l+nreso])**2)/para_s0[i,l+2*nreso]**2)               
            ampm1[k]+= para_s1[i,l]/(1+((para_s1[i,nreso+k]-\
                para_s1[i,l+nreso])**2)/para_s1[i,l+2*nreso]**2)   
    salida0[i,:nreso]=ampm0[:]
    salida1[i,:nreso]=ampm1[:]
                
        
path0 = pathmonth+"mediaLO"+'_0'
path1 = pathmonth+"mediaLO"+'_1'
pathe0 = pathmonth+"mediaLOE"+'_0'
pathe1 = pathmonth+"mediaLOE"+'_1'
pathc0 = pathmonth+"mediaLOC"+'_0'
pathc1 = pathmonth+"mediaLOC"+'_1'
# OUTPUTS
np.savetxt(path0,salida0.ravel())
np.savetxt(path1,salida1.ravel())
np.savetxt(pathe0,salerr0.ravel())
np.savetxt(pathe1,salerr1.ravel())            
np.savetxt(pathc0,salchi0.ravel())
np.savetxt(pathc1,salchi1.ravel())
    

