#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for preparing SWAN results for assimilation 
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
from numpy import ma
import numpy as np
import datetime as datetime
import netcdftime
import pylab as pl
import sys,os

#swan2netcdf
arg    = sys.argv
dir_in = arg[1]
##############################################
try:
    os.system('rm base_info.pyc'  )
except:
    pass
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
##############################################
inp_dir      = local_inp = base_info.base_dir+'/inp'
real_data    = base_info.real_data
hisfile_name = base_info.hisfile_name
base_dir     = base_info.base_dir
##################################################
#which time is the wave assimilation time?
if real_data:
    #obs_file=inp_dir+'/obs_wave/rad/cBathy2012-05-10T05.nc'
    obs_file = inp_dir+'/obs_wave/rad/cBathy.nc'

    time_var = 'time'
else:
    obs_file = inp_dir+'/obs_wave/syn/wav_syn_obs.nc'
    time_var = 'ocean_time'
###
#if real:
nc       = netCDF4.Dataset(obs_file)
ncv_obs  = nc.variables
xr       = ncv_obs['x'][:]
yr       = ncv_obs['y'][:]
utim_obs = netcdftime.utime(ncv_obs[time_var].units)
sec_obs  = ncv_obs[time_var][:]
date_obs = utim_obs.num2date(sec_obs)
##>>>> Choose which obs data going to use for assimilation.

try:
    date_assim = date_obs[0]
    tim        = sec_obs [0]
except:
    date_assim = date_obs
    tim        = sec_obs

units=ncv_obs[time_var].units
nc.close()

#read bathy file size
local_inp = base_info.base_dir+'/inp'
ncf   = local_inp + '/const/' + base_info.grd
#ncf  = base + '/03_mem_inp/member1004/'+hisfile_name
nc    = netCDF4.Dataset(ncf)
ncv   = nc.variables
bathy = ncv['h'][:]
nc.close()
row,col = bathy.shape

outfile = dir_in+'/swout.nc'
print outfile,date_assim

## read swan data
# preallocate variables for speed
xp     = pl.loadtxt(dir_in + '/xp.blkout'    )
yp     = pl.loadtxt(dir_in + '/yp.blkout'    )
hsign  = pl.loadtxt(dir_in + '/swh.blkout'   )
pdir   = pl.loadtxt(dir_in + '/pdir.blkout'  )
tm01   = pl.loadtxt(dir_in + '/tm01.blkout'  )
depth  = pl.loadtxt(dir_in + '/depth.blkout' )
botlev = pl.loadtxt(dir_in + '/bottom.blkout')
force  = pl.loadtxt(dir_in + '/force.blkout' )
dir    = pl.loadtxt(dir_in + '/dir.blkout'   )
rtp    = pl.loadtxt(dir_in + '/rtp.blkout'   )
dissip = pl.loadtxt(dir_in + '/dissip.blkout')
ubot   = pl.loadtxt(dir_in + '/ubot.blkout'  )
wlen   = pl.loadtxt(dir_in + '/wlen.blkout'  )
qb     = pl.loadtxt(dir_in + '/qb.blkout'    )
transp = pl.loadtxt(dir_in + '/wp.blkout'    )
vel    = pl.loadtxt(dir_in + '/vel.blkout'   )


if False:
    figure()
    teta=270-dir
    teta[teta<360]=teta[teta<360]+360
    hsign=ma.masked_where(hsign==-9.0,hsign)
    hsx=hsign*cos(teta)
    hsy=hsign*sin(teta)
    pcolor(xp,yp,dir)
    kk=3
    quiver(xp[::kk,::kk],yp[::kk,::kk],hsx[::kk,::kk],hsy[::kk,::kk],scale=0.8)



# clean up swan data
# clean up the variables
mask = hsign.copy()
mask[mask>-5]=1.0
mask[mask<-5]=0.0


hsign [hsign <-5.]    = 0.0;
pdir  [pdir  <-5.]    = 0.0;
tm01  [tm01  <-5.]    = 0.0;
depth [depth <-5.]    = 0.0;
botlev[botlev<-5.]    = 0.0;
force [force <-5.]    = 0.0;
dir   [dir   <-5.]    = 0.0;
rtp   [rtp   <-5.]    = 0.0;
dissip[dissip<-5.]    = 0.0;
wlen  [wlen  <-5]     = 0.0;
qb    [qb    <-1]     = 0.0;

# decompose vectors
velx    = vel   [0:row ,:];
vely    = vel   [row:  ,:];
forcex  = force [0:row ,:];
forcey  = force [row:  ,:];
transpx = transp[0:row ,:];
transpy = transp[row:  ,:];

hsign =np.ma.masked_where(mask ==0.0, hsign)
pdir  =np.ma.masked_where(mask ==0.0, pdir)
tm01  =np.ma.masked_where(mask ==0.0, tm01)
depth =np.ma.masked_where(mask ==0.0, depth)
botlev=np.ma.masked_where(mask ==0.0, botlev)
dir   =np.ma.masked_where(mask ==0.0, dir)
rtp   =np.ma.masked_where(mask ==0.0, rtp)
dissip=np.ma.masked_where(mask ==0.0, dissip)
wlen  =np.ma.masked_where(mask ==0.0, wlen)
velx  =np.ma.masked_where(mask ==0.0, velx)
vely  =np.ma.masked_where(mask ==0.0, vely)
forcex=np.ma.masked_where(mask ==0.0, forcex)
forcey=np.ma.masked_where(mask ==0.0, forcey)
transpx  =np.ma.masked_where(mask ==0.0, transpx)
transpy  =np.ma.masked_where(mask ==0.0, transpy)

#
if True:
    #Cut to smaller region
    i1,i2,j1,j2=pl.loadtxt('param.inp')
    i1=int(i1)
    i2=int(i2)
    j1=int(j1)
    j2=int(j2)
    print i1,i2,j1,j2
    
    #cut to smaller region
    #i1=60
    #i2=150
    #j1=30
    #j2=180
    k=1
    xp=xp [j1:j2:k,i1:i2:k]
    yp=yp [j1:j2:k,i1:i2:k]
    hsign=hsign [j1:j2:k,i1:i2:k]
    pdir=pdir [j1:j2:k,i1:i2:k]
    tm01=tm01 [j1:j2:k,i1:i2:k]
    depth=depth[j1:j2:k,i1:i2:k]
    botlev=botlev [j1:j2:k,i1:i2:k]
    forcex=forcex [j1:j2:k,i1:i2:k]
    forcey=forcey [j1:j2:k,i1:i2:k]
    dir=dir [j1:j2:k,i1:i2:k]
    rtp=rtp [j1:j2:k,i1:i2:k]
    dissip=dissip[j1:j2:k,i1:i2:k] 
    ubot=ubot [j1:j2:k,i1:i2:k]
    wlen=wlen [j1:j2:k,i1:i2:k]
    qb=qb [j1:j2:k,i1:i2:k]
    transpx=transpx [j1:j2:k,i1:i2:k]
    transpy=transpy [j1:j2:k,i1:i2:k]
    velx=velx [j1:j2:k,i1:i2:k]
    vely=vely [j1:j2:k,i1:i2:k]
    mask=mask [j1:j2:k,i1:i2:k]
    
#print('writing netcdf file')
##________netcdf writing________________

ny,nx = pl.shape(xp)

missing_value=0
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')

nc.createDimension('x', nx)
nc.createDimension('y', ny)

nc.createDimension('ocean_time',None)
timea = nc.createVariable('ocean_time','f8',('ocean_time',))
timea.units = units
timea[:]=tim

x_nc = nc.createVariable('x', 'float', ('y','x',))
x_nc.long_name = 'x positions'
x_nc[:] = xp[:]
 
y_nc = nc.createVariable('y', 'float', ('y','x',))
y_nc.long_name = 'y positions'
y_nc[:] = yp[:]  #[:,1]
    
p10 = nc.createVariable('hsig', 'float', ('ocean_time','y','x',))
p10.long_name = 'significant wave height'
p10.units = 'm'
p10[0,:] = hsign[:]

p11 = nc.createVariable('pdir', 'float', ('ocean_time','y','x',))
p11.long_name = 'peak wave direction'
p11.units = 'deg'
p11[0,:] = pdir[:]

p12 = nc.createVariable('tm01', 'float', ('ocean_time','y','x',))
p12.long_name = 'mean wave period'
p12.units = 's'
p12[0,:] = tm01[:]

p13 = nc.createVariable('dep', 'float', ('ocean_time','y','x',))
p13.long_name = 'water depth'
p13.units = 'm'
p13[0,:] = depth[:]

p14 = nc.createVariable('forcex', 'float', ('ocean_time','y','x',))
p14.long_name = 'wave-induced force x'
p14.units = 'nm-2'
p14[0,:] = forcex[:]

p142 = nc.createVariable('forcey', 'float', ('ocean_time','y','x',))
p142.long_name = 'wave-induced force y'
p142.units = 'nm-2'
p142[0,:] = forcey[:]

p15 = nc.createVariable('dir', 'float', ('ocean_time','y','x',))
p15.long_name = 'mean wave direction'
p15.units = 'deg'
p15[0,:] = dir[:]

p16 = nc.createVariable('rtp', 'float', ('ocean_time','y','x',))
p16.long_name = 'peak wave period'
p16.units = 's'
p16[0,:] = rtp[:]

p17 = nc.createVariable('dissip', 'float', ('ocean_time','y','x',))
p17.long_name = 'dissipation'
p17.units = 'wm-2'
p17[0,:] = dissip[:]

p18 = nc.createVariable('ubot', 'float', ('ocean_time','y','x',))
p18.long_name = 'orbital velocity'
p18.units = 'ms-1'
p18[0,:] = ubot[:]

p19 = nc.createVariable('wlen', 'float', ('ocean_time','y','x',))
p19.long_name = 'wave length'
p19.units = 'm'
p19[0,:] = wlen[:]

p21 = nc.createVariable('qb', 'float', ('ocean_time','y','x',))
p21.long_name = 'franction of breaking waves'
p21.units = '0 ~ 1'
p21[0,:] = qb[:]

p22 = nc.createVariable('velx', 'float', ('ocean_time','y','x',))
p22.long_name = 'current velocity x'
p22.units = 'ms-1'
p22[0,:] = velx[:]

p23 = nc.createVariable('vely', 'float', ('ocean_time','y','x',))
p23.long_name = 'current velocity y'
p23.units = 'ms-1'
p23[0,:] = vely[:]

p24 = nc.createVariable('mask', 'float', ('ocean_time','y','x',))
p24.long_name = 'mask hs'
p24[0,:] = mask[:]

nc.created = datetime.datetime.now().isoformat()
nc.close()

