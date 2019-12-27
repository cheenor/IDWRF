#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:46:55 2019

@author: chenjh
"""
import matplotlib as mpl
#mpl.use('Agg')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from wrf import (to_np, getvar, latlon_coords,get_basemap)
#import time_lon_loop2
import numpy.f2py
import numpy as np
from pylab import *
#####################
def getcomdbz(dbz):
    nz=len(dbz[:,0,0])
    ny=len(dbz[0,:,0])
    nx=len(dbz[0,0,:])
    print 'nx,ny',nx,ny
    cob=np.zeros(shape=(ny,nx),dtype=float)
    for iy in range(0,ny):
        for ix in range(0,nx):
            tmp=dbz[:,iy,ix]
            cob[iy,ix]=max(tmp)
    return cob
#####################################################################
dirin='/glade/u/home/chenjh/workdir/IdealWRF/T/WRFV3/run/'
dirout='/glade/u/home/chenjh/workdir/IdealWRF/T/Pics/'
dirin='/Volumes/DATA02/ModelOUTPUT/IDWRF/'
dirin='/Volumes/MacHD/Work/Papers/ICMW2020/Data/'
dirout='/Volumes/MacHD/Work/Papers/ICMW2020/Pics/'
#
#ny=99
#nx=99
#xdat=np.zeros(shape=(nx),dtype=float)
#ydat=np.zeros(shape=(ny),dtype=float)
#for i in range(0,nx):
#    xdat[i]=i*0.1
#for i in range(0,ny):
#    ydat[i]=i*0.1    
fold='R01'
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 14,
        }
nz=80
for dm in range(0,1):
    dmstr='R'+"%2.2d"%(dm)
    for ih in range(0,4):
        hourstr="%2.2d"%ih
        for imn in range(0,6,1):
            minstr="%2.2d"%(imn*10)
            fname='wrfout_d01_0001-01-01_'+hourstr+':'+minstr+':00'
            print fname
            datestr=hourstr+'+'+minstr       
            fpath=dirin+fold+'/'+fname
            wrf_file = Dataset(fpath)
            imm=-1           
            #qv=   wrf_file.variables['QVAPOR'] #getvar(wrf_file, "QVAPOR",timeidx=imm)
            qv=getvar(wrf_file, "QVAPOR",timeidx=imm)
            qc=   getvar(wrf_file, "QCLOUD",timeidx=imm)
            qr=    getvar(wrf_file, "QRAIN",timeidx=imm)
            zdat=   getvar(wrf_file, "z",timeidx=imm,untis='km')
            ny=len(qr[10,:,0])
            nx=len(qr[10,0,:])
            print 'nx,ny=',nx,ny
            ny1=int(ny/2)-10
            ny2=int(ny/2)+10
            xdat=np.zeros(shape=(nx),dtype=float)
            zdat=np.zeros(shape=(nz),dtype=float)
            qcld=np.zeros(shape=(nz,nx),dtype=float)
            qrain=np.zeros(shape=(nz,nx),dtype=float)
            for iz in range(0,nz):
                for ix in range(0,nx):
                    qcld[iz,ix]=0.0
                    qr[iz,ix]=0.0
                    for iy in range(ny1,ny2):
                        qcld[iz,ix]=qcld[iz,ix]+qc[iz,iy,ix]/(ny2-ny1)
                        qrain[iz,ix]=qr[iz,ix]+qr[iz,iy,ix]/(ny2-ny1)
            fig, axs = plt.subplots(nrows=2,ncols=1, figsize=(10,10))
            plt.subplot(2,1,1)
            axs[0]=plt.contourf(xdat,zdat,qcld,
                cmap=plt.get_cmap('Greens'), extend='both')
            marknm='(a) '+fold+' cloud water mixing ratio'
            plt.title(marknm,fontsize=14)
            plt.title(marknm,fontsize=14)        
            axx=plt.subplot(2,1,1)
            xmajorLocator   = MultipleLocator(2)
            axx.xaxis.set_major_locator(xmajorLocator) 
            ymajorLocator   = MultipleLocator(2)
            axx.yaxis.set_major_locator(ymajorLocator)
            plt.xlabel(r'X ($km$)', fontdict=font)
            plt.ylabel(r'Z ($km$)', fontdict=font)
            plt.subplot(2,1,2)
            axs[1]=plt.contourf(xdat,zdat,qrain,
                cmap=plt.get_cmap('Greens'), extend='both')
            marknm='(a) '+fold+' rain water mixing ratio'
            plt.title(marknm,fontsize=14)        
            axx=plt.subplot(2,1,2)
            xmajorLocator   = MultipleLocator(2)
            axx.xaxis.set_major_locator(xmajorLocator) 
            ymajorLocator   = MultipleLocator(2)
            axx.yaxis.set_major_locator(ymajorLocator)
            plt.xlabel(r'X ($km$)', fontdict=font)
            plt.ylabel(r'Z ($km$)', fontdict=font)
            plt.subplots_adjust(left = 0.08, wspace = 0.2, hspace = 0.4, \
                bottom = 0.05, right=0.9, top = 0.95)
            plt.colorbar(orientation='horizontal')
            plt.savefig(dirout+fold+"_dbz_"+datestr+".png",dpi=300)          
            plt.close()



