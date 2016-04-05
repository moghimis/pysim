#!/usr/bin/python
#==============================================================#
__title__          = "test_generate2D.py"
__description__    = "Script to test generate2D.py"
__author__         = "Cigdem Akan"
__date__           = "22-June-2015"
__email__          = "akanc@ucla.edu"
__pyhton_version__ = "2.7.9"
#--------------------------------------------------------------#
# NOTES:
# 1- Uses method of Evensen (1994) for generating random 2D
#    fields. The output is an ensemble of N realizations,
#    having covariance
#
#       C(dx,dy) = Lz*exp(-3*((dx/Lx)^2+(dy/Ly)^2))
#
#    This means the covariance is roughly zero at distances
#    greater than |Lx,Ly|, and the covariance is equal to
#    exp(-1) at about |Lx,Ly|/3.
#
# 2- Adopted from G. Wilson's Matlab function.
#--------------------------------------------------------------#
# LOGS:
# 1- 19-June-2015: For some reason fsolve is not working properly
#                  when the input is a vector. So, I had to use
#                  p and l instead of pn0, ln0.
#==============================================================#
import os
import sys
import scipy as sp
import scipy.sparse
from numpy import *
import glob
import pylab as pl
import matplotlib as plt
import matplotlib.colors as pltc
import matplotlib.pyplot as pylt
import numpy as np
from scipy import interpolate
import scipy.optimize
import cmath
import numpy.matlib as npm
import scipy.io as sio
import netCDF4
from generate2D import generate2D

# Read the netcdf file
grd   = netCDF4.Dataset('prior.nc')
grdv  = grd.variables
h     = grdv['h'][:]
xrho  = grdv['x_rho'][:]
maskr = grdv['mask_rho'][:]
dx    = grdv['pm'][:]
dx    = 1./np.min(dx)
dy    = grdv['pn'][:]
dy    = 1./np.min(dy)
nx,ny = xrho.shape

# Input parameters
N     = 100
pLx   = 200.
pLy   = 200.
pLz   = 0.5

dhN  = generate2D(nx,ny,dx,dy,pLx,pLy,pLz,N);
dh   = dhN[:,:,N-1]
minh = np.nanmin(h)
maxh = np.nanmax(h)

#h[h<minh] = 0.2
h[maskr==0]=nan

dc  = 65
vz  = np.linspace(-1,12,dc,endpoint=True)
vz0 = np.linspace(0.,0.5,dc,endpoint=True)
 
fig = pylt.gcf()
fig.set_size_inches(12,12)
pylt.subplot(2,2,1)
pylt.contourf(h,vz,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('h')
#
pylt.subplot(2,2,2)
pylt.contourf(h+dh,vz,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('h+dh')
#
pylt.subplot(2,2,3)
pylt.contourf(np.nanmean(dhN,axis=2),cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('mean dhN')
#
pylt.subplot(2,2,4)
pylt.contourf(np.nanstd(dhN,axis=2),cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('std dhN')
pylt.show()

