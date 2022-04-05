#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:02:50 2020
Revised on 2022/04/01
@author: Alfonso
"""
# INPUTS
year='2015'
month='1503'
pathin='S_N_Data/'
pathout='S_N_FD/'
pin=pathin+year+'/'+month+'/'
pou=pathout+year+'/'+month+'/'

import numpy as np
divi=6
nleeventa=2560  # 10 s SIGNAL
nventa=2**13    # df window
nventa2=nventa//2+1
fmues=256       # Hz sampling freq
facdat=10.0/2**15   # Sampling factor
df=float(fmues)/float(nventa)
nhora=fmues*3600  # 921600   # Hourly data
limsup=9.990        # Saturation superior bound
liminf=-9.990        # Saturation lower bound
(f1,f2)=(6,25)   # Bands limits 
nhoradivi=nhora//divi
nw=nhoradivi          # Whel's method
mw=nleeventa
pw=nleeventa//2
ninter=int(np.round((nw-mw)/pw+1))
def hann(m):                          # Hann window
    val=np.arange(0,m)
    hannv=1.0-np.cos(val*np.pi/m)**2
    nor=np.sqrt(np.sum(hannv**2)/m)
    return hannv/nor
hannvs=hann(mw)
datosventa=np.zeros(nventa,dtype=float)     # For fft
datosDF=np.zeros(nventa2,dtype=float)
def pofre(f,df,fi):
    return int(round((f-fi)/df))
listafrec=np.arange(0,nventa2)*fmues/nventa   # frequency list
listafreccal=listafrec[pofre(f1,df,0.):pofre(f2,df,0.)+1]       # calibrated frequec

# Calibration
Lmag= 298E3
Cs= 40E-12
Cc= 0.65*55E-12    #35.75E-12
Ce= 2E-12
Ctot= Cs+Cc+Ce
Rmag= 320E3
facAmp= 2500.0
facSen= 1.9E6
f0= 1/(2.0*np.pi*np.sqrt(Ctot*Lmag))
delta= Rmag/2.0*np.sqrt(Ctot/Lmag)
def Scal(f):
    fac1= 1.E12/(facAmp*facSen)
    fac2= np.sqrt((1-(f/f0)**2)**2+(2.*delta*f/f0)**2)
    fac3= 1/f
    return fac1*fac2*fac3
funcal= Scal(listafreccal)

# Reading files
ficheros0=open(pin+'ficheros0')
sensor0=[]
for line in ficheros0:
    sensor0.append(line[:-1])
ficheros1=open(pin+'ficheros1')
sensor1=[]
for line in ficheros1:
    sensor1.append(line[:-1])
if (len(sensor0) != len(sensor1)):
    print("Alert in files")
nhoras=len(sensor0)
print("Processing {0} hours".format(nhoras))

# Sensors loop
for ise in range(2):
    print("Sensor {0}".format(ise))
    if (ise == 0):
        leesensor=sensor0
    else:
        leesensor=sensor1
# Output files
    satper=open(pou+'SR'+month+'_satper_'+str(ise),'w')
    media=open(pou+'SR'+month+'_media_'+str(ise),'w')

    for ihora in range(nhoras):
        print("Hour {0}".format(ihora))
# Reading hours
        hora= np.array(np.fromfile(pin+leesensor[ihora],\
                           dtype='int16'),dtype='float')*facdat      
## 
        nsath=len(hora[(liminf>hora)]) + len(hora[(limsup<hora)])
        
        print("{0} saturations have been found ".format(nsath))

# Each hour is divided in 6 parts
        for idivi in range(divi):
            horadivi=hora[idivi*nhoradivi:(idivi+1)*nhoradivi]
            nsat=0
            datosDF[:]=0.
# Each 10 s
            for iinter in range(ninter):
                ninf=iinter*pw
                nsup=ninf+mw
                horainter=horadivi[ninf:nsup]
                nsatinter=len(horainter[(liminf>horainter)]) +\
                          len(horainter[(limsup<horainter)])
                if (nsatinter==0):
                    datosventa[0:mw]=horainter*hannvs
                    datosDF += np.abs(np.fft.rfft(datosventa))**2
                else:
                    nsat += 1
# 10 min interval finished
            if (nsat == ninter):
                print("ALERT: ALL INTERVALS ARE SATURATED")
                print('Sensor {0} hour {1}'.format(ise, ihora+1))
            else:
                datosDF = np.sqrt(datosDF*2./((ninter-nsat)*fmues*mw))
                datosDF[0] /= np.sqrt(2.)
                datosDF[-1] /= np.sqrt(2.)
            datos_Bfield= datosDF[pofre(f1,df,0.):pofre(f2,df,0.)+1]*funcal
            print(*datos_Bfield, sep='\n', file=media)
            print(nsat/ninter, file=satper)

# Hours loop finished
    satper.close()
    media.close()
# Sensors loop finished



