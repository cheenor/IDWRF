#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:52:52 2019

@author: chenjh
"""
import matplotlib
import datetime
import numpy as np
import math
import string
import pymeteo.skewt as skewtpy
import matplotlib.pyplot as plt
import pymeteo.thermo as met
def readAscii(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    timestr=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                try:
                    onedim.append(string.atof(strs))
                except:
                    timestr.append(strs)                        
                #if string.atof(strs)>0 :
                #  print strs
    del linesplit
    f.close()
    return onedim
##########
dirin='/Volumes/MacHD/Work/Papers/ICMW2020/ICMW2020case1setup/'
dirout=dirin
dirpic='/Volumes/MacHD/Work/Papers/ICMW2020/Pics/'
fpath=dirin+'alldata.txt'
onedim=readAscii(fpath,1)
nz=410
pres=np.zeros(shape=(nz),dtype=float)
h=np.zeros(shape=(nz),dtype=float)
t=np.zeros(shape=(nz),dtype=float)
pot=np.zeros(shape=(nz),dtype=float)
qv=np.zeros(shape=(nz),dtype=float)
dew=np.zeros(shape=(nz),dtype=float)
rh=np.zeros(shape=(nz),dtype=float)
u=np.zeros(shape=(nz),dtype=float)
v=np.zeros(shape=(nz),dtype=float)
alldata=np.zeros(shape=(7,nz),dtype=float)
na=int(len(onedim)/nz)
k=0
for iz in range(0,nz):
    for iv in range(0,6):
        k=iz*7+iv+1
        alldata[iv,iz]=onedim[k]
k=7*nz
for iz in range(0,nz):
    for iv in range(6,7):
        kk=iz*2+1+k
        alldata[6,iz]=onedim[kk]
###################
pres[:]=alldata[0,:]
h[:]=alldata[1,:]
t[:]=alldata[2,:]
pot[:]=alldata[3,:]
qv[:]=alldata[4,:]
dew[:]=alldata[5,:]
rh[:]=alldata[6,:]
fpath =dirout+'sounding.txt'
fout=open(fpath,'w')
item='%f '%pres[0]
fout.write(item)
item='%f '%pot[0]
fout.write(item)
item='%f '%qv[0]
fout.write(item)
fout.write('\n')
for iz in range(2,nz):
    item='%f '%h[iz]
    fout.write(item)
    item='%f '%pot[iz]
    fout.write(item)
    item='%f '%qv[iz]
    fout.write(item)
    item='%f '%u[iz]
    fout.write(item)
    item='%f '%v[iz]
    fout.write(item)
    fout.write('\n')
fout.close()
fig, axs = plt.subplots(nrows=1,ncols=3, figsize=(10,10))
a=pot
zdat=h/1000
axs[0].plot(a,zdat,c='r')
axs[0].set_title('Th')
a=qv
axs[1].plot(a,zdat,c='g')
axs[1].set_title('qv')
a=rh
axs[2].plot(a,zdat,c='b')
axs[2].set_title('Rh')
plt.savefig(dirpic+"sounding_t_th_qc_profs.png",dpi=300)          
plt.close()
####
fout=dirpic+'Sounding_t_lnp.png'
skewtpy.myplot_chenjh('ICASE', zdat, pot, pres*100, qv*1e-3, u, v, fout,title='Sounding')
cape=met.CAPE(zdat,pres*100,t+273.15,qv*1e-3,1)
print cape[]
fpath='/Volumes/MacHD/Work/Papers/ICMW2020/Data/'+'input_sounding_bk'
onedim=readAscii(fpath,0)
pres=[]
pot=[]
qv=[]
z=[]
u=[]
v=[]
pres.append(onedim[0])
pot.append(onedim[1])
qv.append(onedim[2])
z.append(0)
u.append(0)
v.append(0)
nz=int((len(onedim)-3)/5)
for iz in range(0,nz):
    k=3+iz*5
    z.append(onedim[k])
    k=3+iz*5+1
    pot.append(onedim[k])
    k=3+iz*5+2
    qv.append(onedim[k])
    k=3+iz*5+3
    u.append(onedim[k])
    k=3+iz*5+4
    v.append(onedim[k])
fig, axs = plt.subplots(nrows=1,ncols=3, figsize=(10,10))
a=pot
zdat=z
axs[0].plot(a,zdat,c='r')
axs[0].set_title('Th')
a=qv
axs[1].plot(a,zdat,c='g')
axs[1].set_title('qv')
a=u
axs[2].plot(a,zdat,c='b')
axs[2].set_title('u')
plt.savefig(dirpic+"wrf_sounding_t_th_qc_profs.png",dpi=300)          
plt.close()
#fout=dirpic+'IC_sounding_txt.png'
#skewtpy.myplot_chenjh('ICASE', zdat, pot, pres*100, qv/1000, u, v, fout,title='Sounding')











