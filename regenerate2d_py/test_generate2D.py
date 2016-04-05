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
from generate2D import generate2D


# Input parameters
dx   = 20.
dy   = 20.
Lx   = 4000.  #doamin size
Ly   = 4000.  #domain size
N    = 3
x    = np.linspace(0,Lx,(Lx+dx)/dx)
y    = np.linspace(0,Ly,(Ly+dy)/dy)
xx,yy = np.meshgrid(x,y)
nx   = len(x)
ny   = len(y)

pLx  = 200.
pLy  = 200.
pLz  = 0.5

dhN  = generate2D(nx,ny,dx,dy,pLx,pLy,pLz,N);
dh   = dhN[:,:,N-1]
minh = 0.25
h    = np.zeros((nx,ny))
      
for i in range (0,nx):
    for j in range (0,ny):
        h[i,j] = 10.

h[h<minh] = 0.25

dc = 65
vz = np.linspace(8,12,dc,endpoint=True)

fig = pylt.gcf()
fig.set_size_inches(12,12)
pylt.subplot(2,2,1)
pylt.contourf(h,vz,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('h')
#
pylt.subplot(2,2,2)
pylt.contourf(dh,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('dh')

pylt.subplot(2,2,3)
pylt.contourf(dhN.mean(2) ,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('mean dhN')

pylt.subplot(2,2,4)
pylt.contourf(dhN.std(2) ,cmap='jet')
pylt.colorbar(format='%.4f')
pylt.title('std dhN')
pylt.show()






