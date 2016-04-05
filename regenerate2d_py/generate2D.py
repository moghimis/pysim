#!/usr/bin/python
#==============================================================#
__title__          = "generate2D.py"
__description__    = "Generates random 2D fields"
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
import cmath
import numpy.matlib as npm
import scipy.io as sio

def generate2D(nx,ny,dx,dy,pLx,pLy,pLz,N):
    # to get a nonperiodic ensemble, define extra "ghost" gridpoints
    n1 = np.round(1.2*nx)
    n2 = np.round(1.2*ny)
    
    n1 = n1+np.mod(n1,2)
    n2 = n2+np.mod(n2,2)
    
    # define constants
    pi2    = 2.0*pi
    deltak = pi2**2./((n1*n2)*dx*dy)
    kappa  = pi2/((n1)*dx)
    kappa2 = kappa**2.
    lmbd   = pi2/((n2)*dy)
    lmbd2  = lmbd**2.
    nreal  = N
     
    # rescale decorrelation lengths such that we will get the
    # following form for the covariance as a function of
    # distance delta:
    #     C(delta)=exp(-3*(delta/Lx)^2)
    
    rx = pLx/np.sqrt(3.0)
    ry = pLy/np.sqrt(3.0) 
    
    #------------------------------------------------------------------
    # solve systems for r1,r2,c
    #------------------------------------------------------------------
    # define wavenumber indeces p,l, excluding p==l==0
    p   = np.linspace((-n2/2.+1.),(n2/2.),(n2/2.)-(-n2/2.+1.)+1)
    l   = np.linspace((-n1/2+1),(n1/2),(n1/2)-(-n1/2+1)+1)
    p,l = np.meshgrid(p,l)
    
    # Commented the following lines due to the problem mentioned in LOGS-1
    pp  = np.array(p).flatten()
    ll  = np.array(p).flatten()
    #ind = sp.setdiff1d(np.linspace(0,p.size-1,p.size-1-0+1),sp.where((p==0) & (l==0)))
    ind = sp.setdiff1d(np.linspace(0,p.size-1,p.size-1-0+1),np.r_[sp.where((p==0) & (l==0))])
    ind = ind.astype(int)
    pn0 = pp[ind]
    ln0 = ll[ind]
    
    def ff(ss):
        r1,r2 = ss
        e     = np.exp(-2.0*(kappa2*(ln0**2.)/(r1**2.) + lmbd2*(pn0**2.)/(r2**2.)))
        f     = np.sum(e*(np.cos(kappa*ln0*rx)-np.exp(-1.)))
        g     = np.sum(e*(np.cos(lmbd*pn0*ry)-np.exp(-1.)))
        return (f,g)
    
    r1,r2 = sp.optimize.fsolve(ff,(3.0/rx,3.0/ry))
    
    summ  = np.sum(np.sum(np.exp(-2.0*(kappa2*(l**2.)/(r1**2.)+lmbd2*(p**2.)/(r2**2.)))))
    summ  = summ-1.0
    c     = np.sqrt(1.0/(deltak*summ))
    
    # define aij matrices.  Note rotation is not enabled in this code
    a11   = 1.0/r1**2
    a22   = 1.0/r2**2
    a12   = 0.0*a11
    
    # define wavenumber indeces following matlab ifft2 convention
    l     = np.linspace(0,(n1/2),(n1/2)-0+1)
    p     = np.linspace(0,(n2/2),(n2/2)-0+1)
    p,l   = np.meshgrid(p,l)
    
    # define amplitudes 'C', in 1st quadrant
    e      = np.exp(-( a11*kappa2*(l**2.) + 2.0*a12*kappa*lmbd*l*p + a22*lmbd2*(p**2.) ))
    C      = e*c*np.sqrt(deltak)
    C[0,:] = 0.
    C[:,0] = 0.
    
    # for each wavenumber (p,l) of each sample (j=1..N)
    A = np.zeros((n1,n2,N))
    for nn in range (0,int(nreal)):
        print "Working on ensemble number " + str(nn)
        qhat  = np.zeros((n1,n2))+0j
        qhat2 = np.zeros((n1,n2))+0j
        # 1st quadrant: phase is arbitrary
        phi   = 2.*pi*np.random.random(C.shape)
        phi[:,int(n2)/2] = 0.
        phi[int(n1)/2,:] = 0.
        qhat[0:int(n1)/2+1,0:int(n2)/2+1] = C*np.exp(cmath.sqrt(-1.)*(phi))
        # 3rd quadrant: phase is also arbitrary
        phi2 = 2.*pi*np.random.random(C.shape)
        phi2[:,int(n2)/2] = 0.
        phi2[int(n1)/2,:] = 0.
        qhat2[0:int(n1)/2+1,0:int(n2)/2+1] = C*np.exp(cmath.sqrt(-1.)*(phi2))
        for j in range (int(n1)/2,int(n1)-1):
            for i in range (0,int(n2)/2):
                qhat[j+1,i+1] = np.conj(qhat2[(int(n1)-j)+1,i+1])
        
        qhat[int(n1)/2:int(n1)-2,1]=0.
        # 2nd and 4th quadrants are set by conjugate symmetry
        for i in range (int(n2)/2+1,int(n2)):
            for j in range (0,int(n1)):
                qhat[j,i] = np.conj(qhat[np.mod(int(n1)-j+1,int(n1)),np.mod(int(n2)-i+1,int(n2)+1)])
        
        #print nn
        # Invert the fourier transform to get the sample
        A[:,:,nn] = np.fft.ifft2(qhat)*n1*n2
    
    # cut down to desired size
    A = A[0:nx,0:ny,:]
    
    # correct mean and variance
    AA  = np.array([np.tile(np.mean(A,axis=2), (1,1)) for ii in xrange(int(N))])
    AA  = AA.transpose((1,2,0))
    A   = A-AA
    del AA
    AA  = np.array([np.tile(np.std(A,axis=2), (1,1)) for ii in xrange(int(N))])
    AA  = AA.transpose((1,2,0))
    A   = A/AA*pLz
    del AA
    
    return A
