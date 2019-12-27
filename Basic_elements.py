#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 4 2019

@author: chenjh
"""
import Floop
import matplotlib as mpl
#mpl.use('Agg')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from wrf import (to_np, getvar, latlon_coords,get_basemap)
#import Floop
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
            dbz = getvar(wrf_file, "dbz",timeidx=imm)
            rain=getvar(wrf_file, "RAINNC",timeidx=imm)
            #tk=  getvar(wrf_file, "tk",timeidx=imm)
            #comdbz=np.zeros(shape=(ny,nx),dtype=float)
            #for iy in range(0,ny):
            #    for ix in range(0,nx):
            #        temp=dbz[:,iy,ix]
            #        comdbz[iy,ix]=temp.max()
            comdbz=Floop.getcomdbz(dbz)
            ny=len(qr[10,:,0])
            nx=len(qr[10,0,:])
            print 'nx,ny=',nx,ny
            xdat=np.zeros(shape=(nx),dtype=float)
            ydat=np.zeros(shape=(ny),dtype=float)
            for i in range(0,nx):
                xdat[i]=i*0.1
            for i in range(0,ny):
                ydat[i]=i*0.1 
            fig, axs = plt.subplots(nrows=1,ncols=1, figsize=(10,10))
            '''
            plt.subplot(3,6,ij)
            axs[0,0]=plt.contourf(xdat,ydat,qv,
                cmap=plt.get_cmap('Greens'), extend='both')
            marknm='(a) '+fold+' mixing ratio'
            plt.title(marknm,fontsize=14)    
            axs[0,1]=plt.contourf(xdat,ydat,qc,
                cmap=plt.get_cmap('Greens'), extend='both')
            marknm='(a) '+fold+' mixing ratio'
            plt.title(marknm,fontsize=14) 
            axs[1,0]=plt.contourf(xdat,ydat,qr,
                cmap=plt.get_cmap('Greens'), extend='both')
            marknm='(a) '+fold+' mixing ratio'
            plt.title(marknm,fontsize=14) 
            '''
            aa=comdbz
            plt.contourf(ydat,xdat,aa,
                cmap=plt.get_cmap('Greens'), extend='both')
            '''
            marknm=fold+' dbz at '+datestr
            plt.title(marknm,fontsize=14)
            axx=plt.subplot(1,1,1)
            xmajorLocator   = MultipleLocator(2)
            axx.xaxis.set_major_locator(xmajorLocator) 
            ymajorLocator   = MultipleLocator(2)
            axx.yaxis.set_major_locator(ymajorLocator)
            plt.xlabel(r'X ($km$)', fontdict=font)
            plt.ylabel(r'Y ($km$)', fontdict=font)
            '''
            plt.subplots_adjust(left = 0.08, wspace = 0.2, hspace = 0.4, \
                bottom = 0.05, right=0.9, top = 0.95)
            plt.colorbar(orientation='horizontal')
            plt.savefig(dirout+fold+"_dbz_"+datestr+".png",dpi=300)          
            plt.close()



