#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for adjusting bathy files for ROMS  #######
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

#######################################################
# Logs:
# 1.0 03/25/2013 02:14:41 PM    prepared for general use
# 2.0 changed to use dh instead of total h
#
#

import netCDF4
#from octant.grid import *
from numpy import *
from datetime import datetime
import os,sys
import glob
import pylab as pl
#from mpl_toolkits.basemap import Basemap

#sys.path.append(os.environ["INP_DIR"]+'/scr/py')
#import roms

k=1

h_min=0.1  # in assimilation region
land=-10.0
arg=sys.argv
if len(arg)< 3 :
   print '############################################################################'
   print 'Please try like ... '
   print 'python scr.py  dir_in dir_out '
   print 'good luck!'
   print '############################################################################'
   sys.exit('oops')

dir_in     = arg[1]
dir_out    = arg[2]
ncf_parent = arg[3]     #'nri_grd_2012.nc'   # adapt to ref bathy

i1,i2,j1,j2,hmin = pl.loadtxt('param.inp')
hmin = float(hmin)
#######################################
# read final bathy param
tnc = netCDF4.Dataset(ncf_parent)
ncv = tnc.variables
xc  = ncv['x_rho'][:]
yc  = ncv['y_rho'][:]
h   = ncv['h'][:]
tnc.close
h    = ma.masked_where(h==-10.0,h)
mask = h.mask
h_tru   = h.copy()
h_tru[mask]=0.0
######################################
# read prior bathy param
pnc  = netCDF4.Dataset('../00_prior/prior.nc')
pncv = pnc.variables
hpre = pncv['h'][:]
pnc.close
hpre[hpre==-10.0] = 0
hpre[mask]      = 0

# build adjustCOV for adjusting COV(hh)
#hmin = 1.5
print 'Hard Coded Hmin!!!!!', hmin,' [m]'
adj_cov = zeros((hpre.shape));
adj_cov [hpre > hmin] = 1.0;
ind     = where(hpre < hmin);
adj_cov [ind] = sin( 0.5 * pi * hpre[ind] / hmin ) **2.0
######################################
dir1  = dir_in +'/'+ 'bat*.nc'
flist = glob.glob(dir1)
flist.sort()

for file in flist:    #file=flist[2]
    #file=flist[0]
    ncf_child_in = file
    ncf_child_out= file.replace('01_bat_inp','02_bat_adj')
    comm='cp '+ncf_parent+'  '+ncf_child_out
    os.system(comm)
    
    ################## topo fine read in #############
    tnc2=netCDF4.Dataset(ncf_child_in,'r')
    x =tnc2.variables['x_rho'][:]
    y =tnc2.variables['y_rho'][:]
    hp = hpre + adj_cov * tnc2.variables['h'][:]
    
    # members adjustments
    hp_member = hp.copy()
    hp_member[hp_member<h_min]=h_min    # we had it before in the form of h_min=0.1
    hp_member[mask]=0.0
    
    alfa=pl.zeros_like(hp_member)
    alfa[j1:j2:k,i1:i2:k]=1.0
    #dep_final= h2 * 0.0
    nrow=10
    
    ms1=linspace(0,1,nrow)
    
    ny,nx = hp.shape
    
    #west
    for jm1 in range (ny): 
       alfa[jm1,0:nrow]=ms1
    
    #east
    for im1 in range (ny):  
       alfa[im1,nx-nrow:nx]=linspace(1,0,nrow)
    
    #south
    for im1 in range (nx): 
       alfa[0:nrow,im1]=ms1
    
    #north
    for im1 in range (nx): 
       alfa[ny-nrow:ny,im1]=linspace(1,0,nrow)
    
    
    #imshow(flipud(alfa))
       
    dep_f= alfa * hp_member + (1.0-alfa ) * h_tru
    #dep_f[dep_f < 0.1]= 0.1
    dep_f[mask]=-10.0
    
    nc1=netCDF4.Dataset(ncf_child_out,'r+')
    ncv1=nc1.variables
    ncv1['h'][:]=dep_f[:]
    ncv1['h'].missing_value=land
    ncv1['h'].valid_min = -1.0
    ncv1['h'].valid_max = 50.0
    nc1.close()

#comm='cp  '+ncf_child_out+'  '+new_prior
#os.system(comm)



"""























































ny,nx=shape(h2)

#tnc3=netCDF4.Dataset(ncf_child_out,'r+')

#if True:
off=10
beta2     = zeros_like(h2)
beta2[off:-off,off:-off]=1.0

dep_f=zeros_like(h2)

dep_f=beta2 * h2+ (1.0 - beta2) * h

###### alfa matrix
#if False:
########################## Preperation of final dep combination ######
alfa     = beta2
#dep_final= h2 * 0.0
nrow=10


ms1=linspace(0,1,nrow)

for jm1 in range (off,ny-off): 
   alfa[jm1,off:nrow+off]=ms1

for im1 in range (off,ny-off): 
   alfa[im1,(nx-nrow-off):nx-off]=linspace(1,0,nrow)

#south
for im1 in range (2*off,nx-2*off): 
   alfa[off:nrow+off,im1]=ms1

#north
for im1 in range (2*off,nx-2*off): 
   alfa[(ny-nrow-off):ny-off,im1]=linspace(1,0,nrow)

#imshow(flipud(alfa))
   
dep_f= alfa * h2 + (1.0-alfa ) * h

hf = dep_f.copy()
smoothing=False
if smoothing:
    print '    > smooth the bathymetry and straighten out the edges'
    for i in range(1, h.shape[1]-1):
        for j in range(1, h.shape[0]-1):
            hf[j, i] = 0.96*h[j, i] + 0.01*(h[j+1, i] + 
                                      h[j-1, i] + h[j, i+1] + h[j, i-1])
    
#print '    > modifying min/max depth'
h_min = 0.1
h_max = 50.0

#    hf[isnan(hf)] = h_min
#    hf[hf<h_min]  = h_min
#    hf[hf>h_max]  = h_max
#   
#    grdf.h = hf.copy()
#  
#    ###  Create mask for grid
#    hf=ma.masked_where(hf<=h_min,hf)
#    grdf.mask[hf.mask]=0
#    grdf.mask[~hf.mask]=1

maskc=zeros_like(hf)
maskc[mask]=0
maskc[~mask]=1
hf=hf * maskc + (maskc-1.0) * 10
hf[isnan(hf)] = h_min
hf[hf<h_min]  = h_min
hf[hf>h_max]  = h_max

grdf.h = hf.copy()
grdf.mask=maskc

print '         >>>>>>>>>> ',ncf_child_out
rom_grd_name=ncf_child_out
roms.write_grd(grdf, rom_grd_name, verbose=False)    

"""
