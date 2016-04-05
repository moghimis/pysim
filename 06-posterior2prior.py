from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for preparing new prior based on recent assimilation 
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "moghimis@gmail.com"
###############################################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 03/25/2013 02:14:41 PM   
#
#
#
################################################################

import netCDF4
from numpy import *
import os,sys
import glob
import pylab as pl
import scipy.io as sio

arg = sys.argv
print arg

inp_file_type = arg[1]
#Cut to smaller region
#i1=60
#i2=150
#
#j1=30
#j2=180
smooth_new_prior = False
i1,i2,j1,j2 = pl.loadtxt('param.inp')
i1=int(i1)
i2=int(i2)
j1=int(j1)
j2=int(j2)
print i1,i2,j1,j2

k    = 1
land = -10.0

ncf_parent = 'final_grd.nc'
new_prior  = 'new_prior.nc'
#comm= 'cp  '+ncf_parent+'  '+new_prior
comm = 'cp  '+ncf_parent+'  tmp.nc'
os.system(comm)
################## topo 200 read in #############
tnc  = netCDF4.Dataset(ncf_parent)
ncv  = tnc.variables
xc   = ncv['x_rho'][:]
yc   = ncv['y_rho'][:]
h    = ncv['h'][:]
mcr  = ncv['mask_rho'][:]
tnc.close

h    = ma.masked_where(h==-10.0,h)
mask = h.mask

h2   = h.copy()
h2[mask]=0.0

if inp_file_type !='netcdf':
    flist=glob.glob('*.mat')
    flist.sort()
    matfile=flist[0]
    data=sio.loadmat(matfile) 
    x=data['x'][:]
    y=data['y'][:]
    hp=data['ha1'][:]
    read_me=data['read_me'][:]
else:
    flist   = glob.glob('assimi*.nc')
    flist.sort()
    ncv_asim = netCDF4.Dataset(flist[0]).variables
    x  = ncv_asim['x_rho'][:]
    y  = ncv_asim['y_rho'][:]
    hp = ncv_asim['h_post'][:].mean(2)
    read_me =  netCDF4.Dataset(flist[0]).history
    
hp_all=pl.zeros_like(h)
hp_all[j1:j2:k,i1:i2:k]=hp
hp_all[mask]=0.0

alfa=pl.zeros_like(hp_all)
alfa[j1:j2:k,i1:i2:k]=1.0
#dep_final= h2 * 0.0
nrow=10

ms1=linspace(0,1,nrow)

#west
for jm1 in range (j1,j2): 
   alfa[jm1,i1:nrow+i1]=ms1
   alfa[jm1,i1:nrow+i1]=ms1

#east
for im1 in range (j1,j2):  
   alfa[im1,i2-nrow:i2]=linspace(1,0,nrow)

#south
for im1 in range (i1,i2): 
   alfa[j1:nrow+j1,im1]=ms1

#north
for im1 in range (i1,i2): 
   alfa[j2-nrow:j2,im1]=linspace(1,0,nrow)

#imshow(flipud(alfa))
   

if smooth_new_prior:
    #imshow(flipud(hp_all))
    import octant.csa as csa
    def interp3(x,y,b,xnew,ynew,method):
        xf=x.flatten()
        yf=y.flatten()
        bf=b.flatten()
        interpm=method
        if interpm=='tri':
            xnewf=xnew.flatten()
            ynewf=ynew.flatten()
            print 'tri interp ...'
            from delaunay import  triangulate
            tri=triangulate.Triangulation(xf, yf)
            interp_b=tri.nn_extrapolator(bf)
            bnewf = interp_b(xnewf,ynewf)
            bnew=bnewf.reshape(xnew.shape)
        elif interpm=='csa' :
            print 'csa interp ...'
            import octant.csa as csa
            csa_interp = csa.CSA(xf, yf,bf)
            bnew = csa_interp(xnew,ynew)
        elif interpm=='grd' :
            print 'griddata interp ...'
            bnew=pl.griddata(xf,yf,bf,xnew,ynew)
        return bnew

    ## Here we smooth and remove shallow water areas
    remove_islands = 1
    every   = 2
    hp_all [hp_all < remove_islands] = remove_islands
    hr0 = hp_all[::every,::every] * mcr [::every,::every]
    xr0 = xc     [::every,::every]
    yr0 = yc     [::every,::every]
    hp_all  = interp3(xr0,yr0,hr0,xc,yc,'csa')

    hf22 = hp_all.copy()
    for i in range(1, hf22.shape[1]-1):
        for j in range(1, hf22.shape[0]-1):
            hf22[j, i] = 0.5*hf22[j, i] + 0.125*(hf22[j+1, i] + 
                                      hf22[j-1, i] + hf22[j, i+1] + hf22[j, i-1])    
    hp_all = hf22


## Apply final mask and adjust to boundaries
dep_f= alfa * hp_all + (1.0-alfa ) * h2


if False:
    def smoothing(hf):
        h2=hf.copy()
        print 'smooth the bathymetry and straighten out the edges'
        buf=4
        for i in range(2, hf.shape[1]-2):
            for j in range(2, hf.shape[0]-2):
                #Saeed invented diagonal smoothing :D
                coef=1./9.
                h2[j, i] = coef*(hf[j, i]+
                                 hf[j+1,i+1]+hf[j+2,i+2]+
                                 hf[j+1,i-1]+hf[j+2,i-2]+
                                 hf[j-1,i+1]+hf[j-2,i+2]+
                                 hf[j-1,i-1]+hf[j-2,i-2])
        return h2
    #######################################################
    
    dep0=dep_f.copy()
    buf=5
    #south
    hf=dep_f[j1-buf:nrow+j1+buf,i1-buf:i2+buf]
    dep_f   [j1-buf:nrow+j1+buf,i1-buf:i2+buf]=smoothing(hf)
    
    ##north
    hf=dep_f[j2-nrow-buf:j2,i1-buf:i2+buf]
    dep_f   [j2-nrow-buf:j2,i1-buf:i2+buf]=smoothing(hf)
    #
    ##east
    j2land=94
    hf=dep_f[j1-buf:j2land+buf,i2-nrow-buf:i2+buf]
    dep_f   [j1-buf:j2land+buf,i2-nrow-buf:i2+buf]=smoothing(hf)
    #
    ##west
    j2land=100
    hf=dep_f[j1-buf:j2land+buf,i1-buf:nrow+i1+buf]
    dep_f   [j1-buf:j2land+buf,i1-buf:nrow+i1+buf]=smoothing(hf)
    #

#dep_f[dep_f < 0.1]= 0.1
dep_f=dep_f+(mcr-1)*20
dep_f[dep_f<-5.0]=land

nc1=netCDF4.Dataset('tmp.nc','r+')
ncv1=nc1.variables
ncv1['h'][:]=dep_f[:]
ncv1['h'].missing_value=land
ncv1['h'].valid_min = -1.0
ncv1['h'].valid_max = 50.0

nc1.assim_info = read_me 
nc1.base_dir   = os.getcwd()
nc1.close()

if inp_file_type !='netcdf':
    nc_name  = matfile.replace('.mat','.nc')
else:
    nc_name  = flist[0]

comm='cp  ' + 'tmp.nc  ' + new_prior
os.system(comm)

#to keep track of stuff
os.system('mkdir     '+ nc_name[:-3])
#os.system('mv  *.*   '+ nc_name[:-3])
#move only files and nor dirs
os.system('find . -maxdepth 1 -type f -exec mv {} '+ nc_name[:-3] +' \;')
os.system('cp  '+nc_name[:-3]+'/'+new_prior+' '+new_prior)

