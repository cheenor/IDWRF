#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 09:33:30 2019

@author: chenjh
"""

import matplotlib as mpl
#mpl.use('Agg')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
from wrf import (to_np, getvar, latlon_coords,get_basemap)
from pylab import *
import pymeteo.skewt as skewtpy
dirin='/Volumes/MacHD/Work/Papers/ICMW2020/Data/'
dirout='/Volumes/MacHD/Work/Papers/ICMW2020/Pics/'
filename='wrfinput_d01_IC'
fpath=dirin+filename
wrf_file = Dataset(fpath)
tk=  getvar(wrf_file, "tk")
qv= getvar(wrf_file, "QVAPOR")
z= getvar(wrf_file, "z",units='m')
th= getvar(wrf_file, "th",units='K')
zdat=z[:,0,0]
fig, axs = plt.subplots(nrows=1,ncols=3, figsize=(10,10))
a=tk[:,0,0]
axs[0].plot(a,zdat,c='r')
a=qv[:,0,0]
axs[0].set_title('Tk')
axs[1].plot(a,zdat,c='g')
a=th[:,0,0]
axs[1].set_title('qv')
axs[2].plot(a,zdat,c='b')
axs[2].set_title('Th')
plt.savefig(dirout+"ic_wrfinput_t_th_qc_profs.png",dpi=300)          
plt.close()
pres= getvar(wrf_file, "p",units='Pa')
u= getvar(wrf_file, "ua",units='m s-1')
v= getvar(wrf_file, "va",units='m s-1')
#mcape= getvar(wrf_file, "cape_3d")
#pot=getvar(wrf_file, "th")
p=pres[:,1,1]
q=qv[:,1,1]
u1=u[:,1,1]
v1=v[:,1,1]
pot=th[:,1,1]
nz=len(pres[:,1,1])
dat=np.zeros(shape=(6,nz),dtype=float)
for iz in range(0,nz):
    dat[0,iz]=z[iz,0,0]
    dat[1,iz]=pres[iz,0,0]
    dat[2,iz]=th[iz,0,0]
    dat[3,iz]=qv[iz,0,0]
    dat[4,iz]=u[iz,0,0]
    dat[5,iz]=v[iz,0,0]
fout=dirout+'Souding_IC.png'
#skewtpy.myplot_chenjh('ICASE', dat[0,:],dat[2,:],dat[1,:],dat[3,:],dat[4,:],dat[5,:], fout,title='Sounding_WRF')
skewtpy.plot('ICASE', dat[0,:],dat[2,:],dat[1,:],dat[3,:],dat[4,:],dat[5,:], fout,title='Sounding_WRF')
