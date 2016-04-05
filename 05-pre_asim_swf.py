from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pre-process SWIFT velocity data and members and put them in vector form.
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "moghimis@gmail.com"

#####################################################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 
#
#
#

import os,sys
import glob
import numpy as np
import netCDF4
from datetime import datetime

args = sys.argv
itr = int(args[1])    

#
os.system('rm base_info.pyc'  )
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
##############################################
base_dir  = base_info.base_dir
inp_dir   = base_info.inp_dir
scr_dir   = base_info.scr_dir
prior     = base_info.prior
final_grd = base_info.grd


inp_dir  = base_dir + '/inp/'
run_id   = '/run_'+str(1000+itr) 

class component:
    def __init__(self,num):
        self.x=np.zeros((num),dtype='float')
        self.y=np.zeros((num),dtype='float')
        self.data=np.zeros((num),dtype='float')
        self.s=np.zeros((num),dtype='float')

swf_member_dir = base_dir + run_id + '/04_swf_adj'

#############################
print 'Loading SWF members '
############################

memdir  = swf_member_dir+'/member1*'
dirlist = np.array(glob.glob(memdir))
dirlist.sort()

nmemc = len(dirlist)

#create variables to be read in
# read SIZE from a sample

#ncf   = dirlist[0]+'/nri_his1000.nc'
ncf   = glob.glob(dirlist[0]+'/nri_h*')[0]
nc    = netCDF4.Dataset(ncf)
ncv   = nc.variables
h     = ncv['h'][:]
xc    = ncv['x'][:]
yc    = ncv['y'][:]
u_obs = ncv['obs_ubar'][:]
v_obs = ncv['obs_vbar'][:]

ndata=len(u_obs)
u = component(num=ndata)
v = component(num=ndata)

# Filling the containers
u.x    = xc
u.y    = yc
u.s    = np.ones_like(xc)
u.data = u_obs

v.x    = xc
v.y    = yc
v.s    = np.ones_like(xc)
v.data = v_obs


meas={}
meas['us'] = u
meas['vs'] = v

#### Memebers info 
uc    = np.zeros((ndata,nmemc),dtype='float')
vc    = np.zeros_like(uc)

failind=[]
for i in range(nmemc):
    if np.mod(i,50)==0: print ' > SWF member read in ', (i*1.0/nmemc)*100.0, '%'
 
    ncf   = glob.glob(dirlist[i]+'/nri_h*')[0]
    nc=netCDF4.Dataset(ncf)
    ncv=nc.variables
    uc[:,i]=ncv['ubar'][:]
    vc[:,i]=ncv['vbar'][:]
    nc.close()

meas['us'].model = uc
meas['vs'].model = vc

###################################################
#adding a dummy var to make it similar to wave dictionary
for field in ['us','vs']: #,'z']:
    obs   = meas[field]
    obs.f = np.zeros_like(obs.x)
    meas[field] = obs

###################################################
# remove any nan possibly created during interpolation
for field in meas.keys():
    obs=meas[field]
    failind   = np.array(np.where(np.isnan(obs.model.sum(1)))).squeeze()
    keepind   = np.setdiff1d(range(len(obs.x)),failind)
    obs.x     = obs.x    [keepind]
    obs.y     = obs.y    [keepind]
    obs.s     = obs.s    [keepind]
    obs.f     = obs.f    [keepind]
    obs.data  = obs.data [keepind]
    obs.model = obs.model[keepind,:]

###################################################
# a hack to use only part of data
if False:    
    for field in meas.keys():
        obs=meas[field]
        #keepind = obs.y > 1100
        keepind  = obs.y > 0
        obs.x     = obs.x    [keepind]
        obs.y     = obs.y    [keepind]
        obs.s     = obs.s    [keepind]
        obs.f     = obs.f    [keepind]
        obs.data  = obs.data [keepind]
        obs.model = obs.model[keepind,:]

##################################################
# prepare output
##################################################
for field in ['us','vs']: #,'z']:
    obs        = meas[field]
    namep      = '/swf_'+ field+'.nc'
    out_name   = swf_member_dir + namep
    outnc      = netCDF4.Dataset(out_name,'w',format='NETCDF3_CLASSIC')
    dim_data   = field+'_data_num'
    dim_model  = field+'_model_num'
    num_data,num_model = obs.model.shape
    outnc.createDimension(dim_data  , num_data )
    outnc.createDimension(dim_model , num_model)

    p0 = outnc.createVariable(field+'_x','f8',(dim_data,))
    p0.missing_value = -9999.0
    p0[:]            = obs.x

    p1 = outnc.createVariable(field+'_y','f8',(dim_data,))
    p1.missing_value = -9999.0
    p1[:]            = obs.y

    p2 = outnc.createVariable(field+'_s','f8',(dim_data,))
    p2.missing_value = -9999.0
    p2[:]            = obs.s

    p3 = outnc.createVariable(field+'_data','f8',(dim_data,))
    p3.missing_value = -9999.0
    p3[:]            = obs.data

    p4 = outnc.createVariable(field+'_f','f8',(dim_data,))
    p4.missing_value = -9999.0
    p4[:]            = obs.f

    p5 = outnc.createVariable(field+'_model','f8',(dim_data,dim_model))
    p5.missing_value = -9999.0
    p5[:]            = obs.model
    outnc.history  = 'moghimis@gmail.com  current structure for assimilation '+datetime.now().isoformat()
    #outnc.att=  infile[-12:]
    outnc.close()

args = sys.argv
scr_name = args[0]    
scr_dir1 = os.getcwd()
os.system('cp -fr  ' + scr_name + '    ' + swf_member_dir)
os.system('cp -fr  base_info.py        ' + swf_member_dir)


print 'Members CUR data is ready at > ', out_name
## end
##################################################

