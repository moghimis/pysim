#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for finding and preparing right time steps of ROMS for assimilation
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

import netCDF4
import okean.roms.roms as okr
from numpy import ma
import numpy as np
import datetime as datetime
import netcdftime
import pylab as pl
#import octant.csa as csa
import sys,os
import glob

latlon=False

arg=sys.argv
if len(arg)< 2 :
   print '############################################################################'
   print 'Please try like ... '
   print 'python scr.py  dir_in '
   print 'good luck!'
   print '############################################################################'
   sys.exit('oops')

dir_in    = arg[1]
nobs      = int(arg[2])
curv_flag = arg[3]

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
if base_info.real_data:
    ref_grd  = inp_dir+'/obs/sar/uASAR.nc'
    time_key ='time'
else:
    ref_grd  = inp_dir+'/obs/syn/syn1nri_his.nc'
    time_key ='ocean_time'
###############################################
print ref_grd
#Cut to smaller region
i1,i2,j1,j2=pl.loadtxt('param.inp')
i1=int(i1)
i2=int(i2)
j1=int(j1)
j2=int(j2)

print i1,i2,j1,j2

print '>>>>>>>>>>>>  ', ref_grd
#if real:
nc=netCDF4.Dataset(ref_grd)
ncv_obs=nc.variables
if latlon:
    x_rho='lon_rho'
    y_rho='lat_rho'
else:
    x_rho='x_rho'
    y_rho='y_rho'
#xr=ncv_obs[x_rho][:]
#yr=ncv_obs[y_rho][:]
utim_obs=netcdftime.utime(ncv_obs[time_key].units)
sec_obs=ncv_obs[time_key][:]
date_obs=utim_obs.num2date(sec_obs)

#comment by saeed 27 mar 2013
#mr=ncv_obs['u'][0,:].mask

if curv_flag=='curv':
    print 'Curvi-linear velocity component will be add ...'
    slpd=ncv_obs['slope'][:]

##>>>> Choose which obs data going to use for assimilation.
date_assim=date_obs[nobs] 

#else:
#    grd=okr.Grid(ref_grd)
#    xr=grd.use('x_rho')
#    yr=grd.use('y_rho')
#    mr=grd.use('mask_rho')
#    hr=grd.use('h')
#    #zr=grd.use('zeta')
#    mr=ma.masked_where(mr==0,mr)

dirlist=glob.glob(dir_in+'/m*')
dirlist.sort()

print '>>>>>>>>>>>' , dirlist

for idir in dirlist[:]:
    odir=idir.replace('03_mem_inp','04_mem_adj')
    os.system("mkdir -p " + odir)
    print odir
    #for file in ['/nri_his.nc','/nri_avg.nc']:
    files = glob.glob(dirlist[0]+'/nri_h*')
    for file in ['/'+hisfile_name]:
        infile=idir+file
        outfile=odir+file[:-3]+str(nobs+1000)+'.nc'
        print outfile
        #taking vertical coordinate from original file
        nchis=netCDF4.Dataset(infile)
        ncvar=nchis.variables
        
        timename_sim='ocean_time'
        utime_sim=netcdftime.utime(ncvar[timename_sim].units)
        times=ncvar[timename_sim][:] 
        sdates_sim=utime_sim.num2date(times)
        
        #To find closest time step to assimilation time
        sec_assim =utime_sim.date2num(date_assim)
        diff=np.abs(sec_assim-times)
        ind=np.where(diff==diff.min())
        #ntime=pl.array(ind).item()    # Because of an error I am using the next 
        ntime=pl.array(ind).min()
        print 'ntime > ',ntime,'   date_simul > ',sdates_sim[ind],'   date_assim  > ',date_assim
        
        xmem=ncvar[x_rho][:]
        ymem=ncvar[y_rho][:]       
        
        if run_type_flag != '3D':
            ubar = ncvar['ubar'][ntime,:]
            vbar = ncvar['vbar'][ntime,:]
        else:
            ubar = ncvar['u'   ][ntime,-1,:]
            vbar = ncvar['v'   ][ntime,-1,:]    
        
        zeta=ncvar['zeta'][ntime,:]
        h=ncvar['h'][:]
        maskr=ncvar['mask_rho'][:]
        masku=ncvar['mask_u'][:]
        maskv=ncvar['mask_v'][:]
        
        
        maskr=ma.masked_where(maskr==0,maskr)
        
        e = zeta[:]* maskr
        
        # find difference between original mask and this member mask
        # force member to have the same mask and apply it to other values
        
        #comment by saeed 27 mar 2013
        #ind=np.where((maskr==0) & (mr==1))
        ind=np.where((maskr==0) )
        maskr[ind]=1
        
        h=h*maskr
        
        e[maskr.mask]=0.0
        
        u=pl.zeros_like(e)
        v=pl.zeros_like(e)
        
        u_tmp   = ubar   [:] * masku
        v_tmp   = vbar   [:] * maskv
        
        
        # interpolate to rho points (aprox.)
        u[:,1:-1]=(u_tmp[:,:-1]+u_tmp[:,1:])/2.
        u[:,0]=    u_tmp[:,0]
        u[:,-1]=   u_tmp[:,-1]
        u[maskr.mask]=0.0
        
        
        v[1:-1,:]    = (v_tmp[:-1,:]+v_tmp[1:,:])/2.
        v[0   ,:]    =  v_tmp[0  ,:]
        v[-1  ,:]    =  v_tmp[-1 ,:]
        v[maskr.mask]=0.0
        
        if curv_flag=='curv':
            uvd=np.abs(np.squeeze(u * np.cos(slpd)+v * np.cos(np.pi/2.0 - slpd)))
        
        if True:
            #Cut to smaller region
            #i1=60
            #i2=150
        
            #j1=30
            #j2=180
        
            k=1
            h1=h.data    [j1:j2:k,i1:i2:k]
            e1=e.data    [j1:j2:k,i1:i2:k]
            u1=u.data    [j1:j2:k,i1:i2:k]
            v1 =v.data   [j1:j2:k,i1:i2:k]
            if curv_flag=='curv':
	        uvd1=uvd.data[j1:j2:k,i1:i2:k]
            #mask_rho_ref=mr.data[j1:j2:k,i1:i2:k]
            mask_rho=ncvar['mask_rho'][j1:j2:k,i1:i2:k]
            #mask_rho_uni=maskr.data[j1:j2:k,i1:i2:k]
            x1=xmem[j1:j2:k,i1:i2:k]
            y1=ymem[j1:j2:k,i1:i2:k]
         
            times=times[ntime:ntime+1]
            
        #Rotate to x,y grid
        #u1,v1=okcl.rot2d(u,v,ang,inverse=True)
        
        
        #print('Writing NetCDF file')
        ##________NETCDF writing________________
        
        ny,nx = pl.shape(h1)
        
        missing_value=0
        nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        
        nc.createDimension(x_rho, nx)
        nc.createDimension(y_rho, ny)
        
        nc.createDimension('ocean_time',None)
        timea = nc.createVariable('ocean_time','f8',('ocean_time',))
        timea.units = ncvar[timename_sim].units
        timea[:]=times
        
        x_nc = nc.createVariable(x_rho, 'float', (y_rho,x_rho,))
        x_nc.long_name = 'X Positions'
        x_nc[:] = x1[:]
         
        y_nc = nc.createVariable(y_rho, 'float', (y_rho,x_rho,))
        y_nc.long_name = 'Y Positions'
        y_nc[:] = y1[:]  #[:,1]
            
        bathy_nc = nc.createVariable('h', 'float', (y_rho,x_rho,))
        bathy_nc.long_name = 'depth'
        bathy_nc.units = 'm'
        bathy_nc[:] = h1[:]
        
        
#        m_nc = nc.createVariable('mask_rho_ref', 'float', (y_rho,x_rho,))
#        m_nc.long_name = 'depth'
#        m_nc.units = 'm'
#        m_nc[:] = mask_rho_ref[:]
        
        m2_nc = nc.createVariable('mask', 'float', (y_rho,x_rho,))
        m2_nc.long_name = 'depth'
        m2_nc.units = 'm'
        m2_nc[:] = mask_rho[:]
        
#        m32_nc = nc.createVariable('mask_rho_uni', 'float', (y_rho,x_rho,))
#        m32_nc.long_name = 'depth'
#        m32_nc.units = 'm'
#        m32_nc[:] = mask_rho_uni[:]
        
        u_nc = nc.createVariable('ubar', 'float', ('ocean_time',y_rho,x_rho,))
        u_nc.long_name = 'Eastward u m/s'
        u_nc.units = 'ms-1'
        u_nc[0] = u1[:]
        
        v_nc = nc.createVariable('vbar', 'float', ('ocean_time',y_rho,x_rho,))
        v_nc.long_name = 'Northward v m/s'
        v_nc.units = 'ms-1'
        v_nc[0] = v1[:]
        
        z_nc = nc.createVariable('zeta', 'float', ('ocean_time',y_rho,x_rho,))
        z_nc.long_name = 'Surface elevation m'
        z_nc.units = 'm'
        z_nc[0] = e1[:]

        if curv_flag=='curv':
          uv_nc = nc.createVariable('uv_curv', 'float', ('ocean_time','y_rho','x_rho',))
          uv_nc.long_name = 'uv on cuvi-linear grid  m/s'
          uv_nc.units = 'ms-1'
          uv_nc[0] = uvd1[:]
        #
        #for it in range(nt):
        #    z_nc[it,:] = z_cut[it].T
        #    v_nc[it,:] = v_cut[it].T
        #    u_nc[it,:] = u_cut[it].T
        
        
        nc.Created = datetime.datetime.now().isoformat()
        nc.close()
        
        

        
        
