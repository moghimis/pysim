#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for preparing SWIFT files from ROMS results 
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

#####################################################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 03/25/2013 02:14:41 PM   
#
#
#

import netCDF4
#import okean.roms.roms as okr
from numpy import ma
import numpy as np
import datetime as datetime
import netcdftime
import pylab as pl
#import octant.csa as csa
import sys,os
import glob


arg=sys.argv
if len(arg)< 2 :
   print '############################################################################'
   print 'Please try like ... '
   print 'python scr.py  dir_in '
   print 'good luck!'
   print '############################################################################'
   sys.exit('oops')
dir_in=arg[1]


##############################################
try:
    os.system('rm base_info.pyc'  )
except:
    pass
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
##############################################
run_type_flag = base_info.run_type
inp_dir       = local_inp = base_info.base_dir+'/inp'
hisfile_name  = base_info.hisfile_name
###############################################
ref_obs=inp_dir +'/obs_swift/swift_obs.nc'
print arg,ref_obs


nc=netCDF4.Dataset(ref_obs)
ncv_obs=nc.variables
xobs=ncv_obs['x'][:]
yobs=ncv_obs['y'][:]
uobs=ncv_obs['u'][:]
vobs=ncv_obs['v'][:]
hobs=ncv_obs['h'][:]

utim_obs=netcdftime.utime(ncv_obs['time'].units)
sec_obs=ncv_obs['time'][:]
date_obs=utim_obs.num2date(sec_obs)

##>>>>>>
dirlist=glob.glob(dir_in+'/m*')
dirlist.sort()

for idir in dirlist[:]:
    odir=idir.replace('03_mem_inp','04_swf_adj')
    os.system("mkdir -p " + odir)
    print idir
    
    for file in ['/'+hisfile_name]:#,'/nri_avg.nc']:
        infile=idir+file
        outfile=odir+file[:-3]+'.nc'
        #taking vertical coordinate from original file
        nchis=netCDF4.Dataset(infile)
        ncvar=nchis.variables
        
        timename_sim='ocean_time'
        utime_sim=netcdftime.utime(ncvar[timename_sim].units)
        sec_sim=ncvar[timename_sim][:] 
        sdates_sim=utime_sim.num2date(sec_sim)
        
        #To find closest time step to assimilation time
        #the same time ref. for both obs and sim
        sec_obs =utime_sim.date2num(date_obs)
        usim=[]
        vsim=[]
        xsim=[]
        ysim=[]
        tsim=[]
        zsim=[]
        hsim=[]
        
        tobs1=[]
        xobs1=[]
        yobs1=[]
        uobs1=[]
        vobs1=[]
        hobs1=[]

        for it in range(len(sec_obs)):
            obs_step=sec_obs[it]
            diff=np.abs(obs_step-sec_sim)
            ind=np.where(diff==diff.min())
            ntime=pl.array(ind).item()
            #print 'ntime > ',ntime            
            obs_cur_time=utime_sim.num2date(obs_step).isoformat()
            sim_cur_time=sdates_sim[ntime].isoformat()
            if diff.min()>3600: 
                print '>>> skip >  sim:',   sim_cur_time, '      obs:', obs_cur_time
                continue

            if run_type_flag != '3D':
                ubar = ncvar['ubar'][ntime,:]
                vbar = ncvar['vbar'][ntime,:]
            else:
                ubar = ncvar['u'   ][ntime,-1,:]
                vbar = ncvar['v'   ][ntime,-1,:]                
            
            zeta = ncvar['zeta'][ntime,:]
            mu   = ncvar['mask_u']
            mv   = ncvar['mask_v']
            mr   = ncvar['mask_rho']
            xr   = ncvar['x_rho'][:]
            yr   = ncvar['y_rho'][:] 
            hr   = ncvar['h'][:] 
            if (xobs[it]>= xr.max() or xobs[it]<=xr.min() or \
                yobs[it]>= yr.max() or yobs[it]<=yr.min()):
                print '*** skip >  obs exceed sim xy   obs:xy> ',xobs[it],yobs[it]
                continue
            
            uvobs=np.sqrt(uobs[it]*uobs[it]+vobs[it]*vobs[it])
            if (uvobs > 3):
                print '=== skip >  Large velocity in data  u: ',   uobs[it],'  v: ',vobs[it],'  date: ',obs_cur_time
                continue
            
            
            dist2=pl.sqrt ((yr-yobs[it])**2+(xr-xobs[it])**2)
            distm2=dist2.min()
            idx2=pl.where(dist2==distm2)    
            [j],[i]=idx2
            
            if i > 198 or i<2 or j > 198 or j<2:
                print '---- skip >  close to BOU'
                continue
            
            # interpolate stagger grid from u-point and v-point to center of the cell
            uij   = ubar[j,i  ] * mu[j,i  ]
            uim1_j= ubar[j,i-1] * mu[j,i-1]
            
            vij   = vbar[j,i  ] * mv[j,i  ]
            vi_jm1= vbar[j-1,i] * mv[j-1,i]
            
            zz    = zeta[j,i]   * mr[j,i  ]
            hh    = hr  [j,i]   * mr[j,i  ]
            uu= (uij+uim1_j)/2.0
            vv= (vij+vi_jm1)/2.0
            
            xsim.append(xr  [j,i] )
            ysim.append(yr  [j,i])
            tsim.append(sec_sim[ntime])
            usim.append(uu)
            vsim.append(vv)
            zsim.append(zz)
            hsim.append(hh)
            
            tobs1.append(obs_step)
            xobs1.append(xobs[it])
            yobs1.append(yobs[it])
            uobs1.append(uobs[it])
            vobs1.append(vobs[it])
            hobs1.append(hobs[it])

        usim=np.array(usim)
        vsim=np.array(vsim)
        zsim=np.array(zsim)
        hsim=np.array(hsim)
        xsim=np.array(xsim)
        ysim=np.array(ysim)
        tsim=np.array(tsim)
        
        tobs1=np.array(tobs1)
        xobs1=np.array(xobs1)
        yobs1=np.array(yobs1)
        uobs1=np.array(uobs1)
        vobs1=np.array(vobs1)
        hobs1=np.array(hobs1)      
    

        print 'xobs.shape()  > ',xobs1.shape

        
                
        #print('Writing NetCDF file')
        ##________NETCDF writing________________
       
        nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        
        nc.createDimension('ocean_time',None)
        timea = nc.createVariable('ocean_time','f8',('ocean_time',))
        timea.units = ncvar[timename_sim].units
        timea[:]=tsim[:]
        
        x_nc = nc.createVariable   ('x', 'float', ('ocean_time',))
        x_nc.long_name = 'X Positions'
        x_nc[:] = xsim[:]
         
        y_nc = nc.createVariable   ('y', 'float', ('ocean_time',))
        y_nc.long_name = 'Y Positions'
        y_nc[:] = ysim[:]  #[:,1]
            
        bathy_nc = nc.createVariable('h', 'float',('ocean_time',))
        bathy_nc.long_name = 'depth'
        bathy_nc.units = 'm'
        bathy_nc[:] = hsim[:]
       
        u_nc = nc.createVariable('ubar', 'float', ('ocean_time',))
        u_nc.long_name = 'Eastward u m/s'
        u_nc.units = 'ms-1'
        u_nc[:] = usim[:]
        
        v_nc = nc.createVariable('vbar', 'float', ('ocean_time',))
        v_nc.long_name = 'Northward v m/s'
        v_nc.units = 'ms-1'
        v_nc[:] = vsim[:]
        
        z_nc = nc.createVariable('zeta', 'float', ('ocean_time',))
        z_nc.long_name = 'Surface elevation m'
        z_nc.units = 'm'
        z_nc[:] = zsim[:]
        

        ### OBS Part
        timea1 = nc.createVariable('obs_time','f8',('ocean_time',))
        timea1.units = ncvar[timename_sim].units
        timea1[:]=tobs1[:]
        
        x_nc1 = nc.createVariable   ('obs_x', 'float', ('ocean_time',))
        x_nc1.long_name = 'X Positions'
        x_nc1[:] = xobs1[:]
         
        y_nc1 = nc.createVariable   ('obs_y', 'float', ('ocean_time',))
        y_nc1.long_name = 'Y Positions'
        y_nc1[:] = yobs1[:]  #[:,1]
            
        bathy_nc1 = nc.createVariable('obs_h', 'float',('ocean_time',))
        bathy_nc1.long_name = 'total depth m'
        bathy_nc1.units = 'm'
        bathy_nc1[:] = hobs1[:]
       
        u_nc1 = nc.createVariable('obs_ubar', 'float', ('ocean_time',))
        u_nc1.long_name = 'Eastward u m/s'
        u_nc1.units = 'ms-1'
        u_nc1[:] = uobs1[:]
        
        v_nc1 = nc.createVariable('obs_vbar', 'float', ('ocean_time',))
        v_nc1.long_name = 'Northward v m/s'
        v_nc1.units = 'ms-1'
        v_nc1[:] = vobs1[:]
        
        nc.Created = datetime.datetime.now().isoformat()
        nc.history= ' By SM  moghimis@gmail.com'
        nc.close()
        
        

        
        
