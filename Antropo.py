#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 11:39:24 2022

@author: Alfonso
V:2022/04/01
"""
import numpy as np

# INPUTS
year="2015"
month="1503"
pathg='S_N_FD/'           # Have to be specified        

nwindow=2**13
fm=256
df=float(fm)/float(nwindow)
bandcal=np.array([6,25],dtype=float)
nhour=fm*3600
fre=np.arange(bandcal[0],bandcal[1]+df,df)   # calibrated frequec
nf=len(fre)
def pofre(f,df,fi):
    return int(round((f-fi)/df))

pathmonth=pathg+year+"/"+month+"/"+"SR"+month[0:4]+"_"

fichmedia0=pathmonth+"media"+"_0"
fichmedia1=pathmonth+"media"+"_1"

fichsatper0=pathmonth+"satper"+"_0"
fichsatper1=pathmonth+"satper"+"_1"                               
                               
medidas0=np.genfromtxt(fichmedia0).reshape(-1,nf)
medidas1=np.genfromtxt(fichmedia1).reshape(-1,nf)
satper0=np.genfromtxt(fichsatper0)
satper1=np.genfromtxt(fichsatper1)
nintervm=len(medidas0)

def recta(val,frec):
    a=(val[-1]-val[0])/(frec[-1]-frec[0])
    b=val[0]
    valmod=a*(frec-frec[0])+b
    de=np.max((val-valmod)/valmod)+np.min((val-valmod)/valmod)
    return valmod, de

def antro(medidasantro,fre,refpos,antropos,satper,facfil):
    medidas=np.copy(medidasantro)
    for i in range(len(medidas)):
        if satper[i]<0.9999:
            for j in range(len(antropos)):
                desvia = recta(medidas[i,refpos[j,0]:refpos[j,1]+1],\
                               fre[refpos[j,0]:refpos[j,1]+1])[1]
                resul = recta(medidas[i,antropos[j,0]:antropos[j,1]+1],\
                              fre[antropos[j,0]:antropos[j,1]+1])
                if(resul[1]>facfil*np.abs(desvia)):
                    medidas[i,antropos[j,0]:antropos[j,1]+1]=resul[0]
        else:
            print("Saturated interval ",i)  
    return medidas


listafrecantro=np.array([[[14.55,15.35],[13.5,13.9]],\
                         [[16.50,16.90],[15.8,16.4]]])
positions=np.array([pofre(f,df,bandcal[0])\
                   for f in listafrecantro.ravel()],dtype=int)
antropos=(positions.reshape(listafrecantro.shape))[:,0]
refepos=(positions.reshape(listafrecantro.shape))[:,1]

medcor0 = antro(medidas0,fre,refepos,antropos,satper0,2.5)
medcor1 = antro(medidas1,fre,refepos,antropos,satper1,2.5)
fichsal0=open(pathmonth+"mediaNA"+"_0",mode='w')
fichsal1=open(pathmonth+"mediaNA"+"_1",mode='w')

for i in range(nintervm):
    print(*medcor0[i],sep='\n',file=fichsal0)
    print(*medcor1[i],sep='\n',file=fichsal1)

fichsal0.close()
fichsal1.close()
