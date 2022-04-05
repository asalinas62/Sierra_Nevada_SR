#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 19:41:02 2017

@author: Alfonso

V:2022/04/01
"""

pathg="S_N_FD/"
pathdatos='S_N_Data/'
pathgP="S_N_npz/"

import numpy as np
import os
import datetime
def pofre(f,f0,df):
    return int(round((f-f0)/df))

#meses=['1303','1304','1305','1306','1307','1308','1309','1310','1311','1312',
#       '1401','1402','1403','1404','1405','1406','1408','1409','1410','1412',
#       '1501','1502','1503','1504','1505','1506','1507','1508','1509','1510','1511','1512',
#       '1601','1602','1603','1604','1605','1606','1608','1609','1610',
#       '1701','1702']
months=['1503']
divi=6
fajus_inf=6.35  # Defines the lower limit of the fitting band
fajus_sup=23.75 # Defines the upper limit of the fitting band
fcal_inf=6.
fcal_sup= 25.
fre1= 10.0
fre2= 17.0
nventa=2**13
fm=256
df=float(fm)/float(nventa)
nreso=3
fajus_inf_pos= pofre(fajus_inf,fcal_inf,df)
fajus_sup_pos= pofre(fajus_sup,fcal_inf,df)
fre= np.arange(fcal_inf,fcal_sup+df,df)
freajus=fre[fajus_inf_pos:fajus_sup_pos+1]
nfajus=len(freajus)
fre1p=pofre(fre1,fajus_inf,df)
fre2p=pofre(fre2,fajus_inf,df)

def lorentzF(x,par):            
    val=np.zeros(nreso+1,dtype=float)
    for i in np.arange(nreso):
        val[i]=par[i]/(((x-par[nreso+i])/par[nreso*2+i])**2+1)
        val[nreso]=par[nreso*3]*x+par[nreso*3+1]
    return val.sum()

def maxrela(lista):
    ind=0
    x1=lista[0]
    x2=lista[1]
    for i in np.arange(2,len(lista)):
        if x1<x2 and x2>lista[i]:
            ind=i-1
        x1=lista[i-1]
        x2=lista[i]
        return ind

for mes in months:
    anua='20'+mes[:2]
    pathalmes=pathg+"20"+mes[:2]+"/"+mes+"/"
    pathmes=pathalmes+"SR"+mes    
    for sensor in [0,1]:            
        path_sat=pathmes+'_satper_'+str(sensor)
        path_rs=pathmes+'_mediaNA_'+str(sensor)
        path_parLO=pathmes+'_mediaLO_'+str(sensor)
        path_parLOC=pathmes+'_mediaLOC_'+str(sensor)
        path_parLOE=pathmes+'_mediaLOE_'+str(sensor)                
        sat=np.genfromtxt(path_sat)
        nintervmes=len(sat)
        rs=np.genfromtxt(path_rs).reshape(nintervmes,-1)
        rs_ajus=rs[:,fajus_inf_pos:fajus_sup_pos+1]
        parLO=np.genfromtxt(path_parLO).reshape(nintervmes,-1)
        parLOC=np.genfromtxt(path_parLOC)
        parLOE=np.genfromtxt(path_parLOE).reshape(nintervmes,-1)

        rs_LO=np.zeros((nintervmes,nfajus),dtype=float)
        for j in np.arange(nintervmes):
            for i in np.arange(nfajus):
                rs_LO[j,i]=lorentzF(freajus[i],parLO[j,3:])  
                
        pr=np.zeros((nintervmes,nreso*2),dtype=float)
        for j in np.arange(nintervmes):
            bm1= rs_LO[j,0:fre1p]
            pm1= np.argmax(bm1)
            if pm1 == 0 or pm1 == len(bm1)-1 :
                pm1=maxrela(bm1)
            bm2= rs_LO[j,fre1p:fre2p]
            pm2=np.argmax(bm2)
            if pm2 == 0 or pm2 == len(bm2)-1 :
                pm2=maxrela(bm2)
            pm2=pm2+fre1p
            bm3= rs_LO[j,fre2p:]
            pm3= np.argmax(bm3)
            if pm3 == 0 or pm3 == len(bm3)-1 :
                pm3=maxrela(bm3)
            pm3=pm3+fre2p
            pr[j]=np.array([rs_LO[j,pm1],(freajus[pm1]+freajus[pm1-1])/2+\
                            df*(rs_LO[j,pm1]-rs_LO[j,pm1-1])/(2*rs_LO[j,pm1]-\
                            rs_LO[j,pm1-1]-rs_LO[j,pm1+1]),\
                        rs_LO[j,pm2],(freajus[pm2]+freajus[pm2-1])/2+\
                            df*(rs_LO[j,pm2]-rs_LO[j,pm2-1])/(2*rs_LO[j,pm2]-\
                            rs_LO[j,pm2-1]-rs_LO[j,pm2+1]),\
                        rs_LO[j,pm3],(freajus[pm3]+freajus[pm3-1])/2+\
                            df*(rs_LO[j,pm3]-rs_LO[j,pm3-1])/(2*rs_LO[j,pm3]-\
                            rs_LO[j,pm3-1]-rs_LO[j,pm3+1])])       
        
        # TIME
        dt=3.906E+3 # microseconds       
        pathmesdatos=pathdatos+"20"+mes[:2]+'/'+mes+'/'
        fichmesL=sorted(os.listdir(pathmesdatos))
        fichmes=fichmesL.copy()
        for i in fichmesL:
            if i[0]!= 's': 
                fichmes.remove(i)
        nfich=len(fichmes)
        nfich2=(nfich)//2
        if fichmes[0][0] != 's':
            print('Warning, check the data files')
        fichs0dat=fichmes[0:nfich2:2]
        fichs0txt=fichmes[1:nfich2:2]
        fichs1dat=fichmes[nfich2::2]
        fichs1txt=fichmes[nfich2+1::2]
        if sensor == 0:
            ft=fichs0txt
        elif sensor ==1:
            ft=fichs1txt
        else:
            print('Warning, check the sensor identification')
        tiempo=[]
        for fich in ft:
            f=open(pathmesdatos+fich)
            da_ho=f.readlines()[4][22:45]
            day=int(da_ho[0:2])
            month=int(da_ho[3:5])
            year=int(da_ho[6:10])
            hour=int(da_ho[11:13])
            minute=int(da_ho[14:16])
            second=int(da_ho[17:19])
            micro=int(da_ho[20:23]+'000')
            tiho=datetime.datetime(year,month,day,hour,minute,second,micro)
            for i in range(divi):
                intedi=(3600/divi*fm)*dt*i
                timu=tiho+datetime.timedelta(microseconds=intedi)
                tiempo.append(str(timu))
        if(len(tiempo) != nintervmes):
            print('Warning, check the files of the month')
        tiempo=np.array(tiempo, dtype=str)    
                
        pa= 'SN_'+mes+'_'+str(sensor)
        sal= (sat,fre,freajus,rs,rs_LO,parLO,parLOC,parLOE,pr,tiempo)
        np.savez(pathgP+pa,*sal)

