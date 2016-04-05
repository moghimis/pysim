from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Helper functions for calculating dispersion relation with dopller
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"
#
import numpy as np
from   scipy.optimize import fsolve
from scipy.constants import g
#Functions
def approxDispersion(sigma,h):
    # k=approxDispersion(sigma,h)
    # Approximation due to Guo (2002), see discussion here:
    # http://www.johndfenton.com/Papers/Dispersion-Relation.pdf
    # Inputs should be of like dimension
    # h : total dep
    # sigma = 2*pi*f
    k = sigma**2.0/g*(1.0-np.exp(-(sigma*np.sqrt(h/g))**(5.0/2.0)))**(-2.0/5.0);
    return k
#
def exactDispersion(sigma,dep,u,v,ang,kguess):
    """
    Inputs:
     > sigma : 2 * pi * f
     > ang   : relative angle in radians
     > g     : gravity acceleration
     > k_fin : wave number
     >>>     Check the nans or negative wavenumbers on outputs
    """
    ### Solve dispersion relation to calculate wave number
    sol    = lambda kk : (sigma-np.sqrt(u**2.0+v**2.0)*kk*np.cos(ang))**2.0 - g*kk*np.tanh(kk*dep)  
    k_fin = fsolve(sol,kguess) 
    return k_fin
#
def exactDispersionKirby(sigma,dep,u,v,ang,kguess,Hs):
    """
    Inputs:
     > sigma : 2 * pi * f
     > ang   : relative angle in radians
     > Hs    : Wave height
     > g     : gravity acceleration
     > k_fin : wave number
     >>>     Check the nans or negative wavenumbers on outputs
    """
    ### Solve dispersion relation to calculate wave number based on Kirby
    DD = lambda kk,hh : (8.0+np.cosh(4.0*kk*hh)-2.0*np.tanh(kk*hh)**2.0) /(8.0*np.sinh(kk*hh)**4.0)
    f1 = lambda kk,hh : np.tanh(kk*hh)**5.0
    f2 = lambda kk,hh : (kk*hh/np.sinh(kk*hh))**4.0
    ee = lambda kk,HH : kk*HH/2.0;
    Uk = lambda uu,vv,kk,gam : np.sqrt(uu**2.0+vv**2.0)*kk*np.cos(gam);
    #
    sol = lambda k : (sigma-Uk(u,v,k,ang))**2.0 - \
                         g*k*np.tanh( k*dep+f2(k,dep)*ee(k,Hs) ) * \
                         (1.0+f1(k,dep)*ee(k,Hs)**2.0*DD(k,dep))            
    #print sol(0.11)
    
    k_fin=fsolve(sol,kguess)
    return k_fin 
#
def Dispersion(sigma,dep,u,v,ang,kguess,Hs=0.0):
    if Hs.sum() != 0.0:
        k_fin = exactDispersionKirby(sigma=sigma,dep=dep,u=u,v=v,ang=ang,kguess=kguess,Hs=Hs)
    else:                         
        k_fin = exactDispersion     (sigma=sigma,dep=dep,u=u,v=v,ang=ang,kguess=kguess)
    return k_fin



def calc_k(chunks):
    #expect Hs at least as an array of zero
    sigma2,kguess,hi,ui,vi,gamma,Hs=chunks
    k_all=np.zeros_like(kguess) + kguess
    nj,ni=kguess.shape
    
    for j in range(nj):
        #for i in range(ni):
        k_all[j,:] = Dispersion( sigma = sigma2[j,:],       \
                                 dep   = hi[j,:] ,          \
                                 u     = ui[j,:] ,          \
                                 v     = vi[j,:] ,          \
                                 ang   = gamma[j,:] ,       \
                                kguess = kguess[j,:],       \
                                Hs     = Hs[j,:])
    k_all=np.array(k_all)
    return  k_all


def exactDispersionKirby_solve4dep(sigma,k,u,v,ang,dep_guess,Hs):
    """
    Inputs:
     > sigma : 2 * pi * f
     > ang   : relative angle in radians
     > Hs    : Wave height
     > g     : gravity acceleration
     > k_fin : wave number
     >>>     Check the nans or negative wavenumbers on outputs
    """
    ### Solve dispersion relation to calculate wave number based on Kirby
    DD = lambda kk,hh : (8.0+np.cosh(4.0*kk*hh)-2.0*np.tanh(kk*hh)**2.0) /(8.0*np.sinh(kk*hh)**4.0)
    f1 = lambda kk,hh : np.tanh(kk*hh)**5.0
    f2 = lambda kk,hh : (kk*hh/np.sinh(kk*hh))**4.0
    ee = lambda kk,HH : kk*HH/2.0;
    Uk = lambda uu,vv,kk,gam : np.sqrt(uu**2.0+vv**2.0)*kk*np.cos(gam);
    #
    sol = lambda dep : (sigma-Uk(u,v,k,ang))**2.0 - \
                         g*k*np.tanh( k*dep+f2(k,dep)*ee(k,Hs) ) * \
                         (1.0+f1(k,dep)*ee(k,Hs)**2.0*DD(k,dep))            
    #print sol(0.11)
    
    dep_fin = fsolve(sol,dep_guess)
    return dep_fin 

def exactDispersion_solve4dep(sigma,k,u,v,ang,dep_guess):
    """
    Inputs:
     > sigma : 2 * pi * f
     > ang   : relative angle in radians
     > g     : gravity acceleration
     > k_fin : wave number
     >>>     Check the nans or negative wavenumbers on outputs
    """
    ### Solve dispersion relation to calculate wave number
    sol    = lambda dep : (sigma-np.sqrt(u**2.0+v**2.0)*k*np.cos(ang))**2.0 - g*k*np.tanh(k*dep)  
    dep_fin = fsolve(sol,dep_guess) 
    return dep_fin


if __name__ == "__main__":

    #data    
    sigma0  = 2 * np.pi / 6.
    dep0    = 4.0
    u0      = 0.0
    v0      = 0.0
    ang0    = 0.0
    kguess0 = 0.1
    Hs0     = 1.2
    dep_g0  = 1.0
    
    
    k0      = exactDispersionKirby(sigma=sigma0,dep=dep0,u=u0,v=v0,ang=ang0,kguess=kguess0,Hs=Hs0)
    print 'k0 from kirby', k0
    
    
    
    k_l = exactDispersion     (sigma=sigma0,dep=dep0,u=u0,v=v0,ang=ang0,kguess=kguess0)
    print k_l
    #
    #
    dep_k = exactDispersionKirby_solve4dep(sigma=sigma0,k=k0,u=u0,v=v0,ang=ang0,dep_guess=dep_g0,Hs=Hs0)
    dep_l = exactDispersion_solve4dep     (sigma=sigma0,k=k0,u=u0,v=v0,ang=ang0,dep_guess=dep_g0)

    print "kirby  depth", dep_k
    print "linear depth", dep_l
    
    np.abs(dep_k - dep_l) / dep_k


    