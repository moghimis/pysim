#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script for prepering SWAN input files from ROMS results  #######
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

##############################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 03/25/2013 02:14:41 PM   
#
#
#
##############################################
import netCDF4
import netcdftime
import pylab as pl
import numpy as np
from datetime import datetime
import os,sys
import glob
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
final_grd    = local_inp+'/const/'+base_info.grd
hisfile_name = base_info.hisfile_name
##################################################
################# Funcs ##########################
#write scaler field in file to be read by swan by idla=1
def wrt_field(outfile,name):
    outf=open(outfile,'w')
    ny,nx=name.shape
    aa=pl.arange(ny-1,-1,-1)
    for iy in aa:
        str1=' '.join('%.4f' % n for n in name[iy,:])
        str1=str1+'\n'
        outf.write(str1)
    outf.close()

#write vector field in file to be read by swan by idla=1
def wrt_vec(outfile,name1,name2):
    outf=open(outfile,'w')
    ny,nx=name1.shape
    aa=pl.arange(ny-1,-1,-1)
    for iy in aa:
        str1=' '.join('%.4f' % n for n in name1[iy,:])
        str1=str1+'\n'
        outf.write(str1)
    
    for iy in aa:
        str1=' '.join('%.4f' % n for n in name2[iy,:])
        str1=str1+'\n'
        outf.write(str1)    
    outf.close()
######################################################
#N=100

if real_data:
    print 'real data assimilation is considered'
    #obs_file=inp_dir+'/obs_wave/rad/cBathy2012-05-10T05.nc'
    obs_file=inp_dir+'/obs_wave/rad/cBathy.nc'

    time_var='time'
else:
    print 'Syn data assimilation is considered'
    obs_file=inp_dir+'/obs_wave/syn/wav_syn_obs.nc'
    time_var='ocean_time'

#if real:
nc       = netCDF4.Dataset(obs_file)
ncv_obs  = nc.variables
xr       = ncv_obs['x'][:]
yr       = ncv_obs['y'][:]
utim_obs = netcdftime.utime(ncv_obs[time_var].units)
sec_obs  = ncv_obs[time_var][:]
date_obs = utim_obs.num2date(sec_obs)
nc.close()
##>>>> Choose which obs data going to use for assimilation.
try:
    date_assim = date_obs[0]
except:
    date_assim = date_obs
  
#ntime=12  # netcdf index
arg    = sys.argv
dirin = arg[1]

#dirin ='../03_mem_inp/member'+ str(1000+n+1)
dirout = dirin.split('/')[-1]+'/'  
print dirout
## Read the ROMS Data
# This is a code to manipulate the input from the NetCDF files given by
# Saeed's ROMS model.
if True:
    # Enter path to the data file
    #datapath = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/seq1/every2/run_1000/03_mem_inp/member1001';
    datapath = dirin;
    datafile = datapath+'/'+hisfile_name;
    tnc      = netCDF4.Dataset(datafile)
    x        = tnc.variables['x_rho'][:]
    y        = tnc.variables['y_rho'][:]
    h        = tnc.variables['h'][:]

    ####   >  find ntme based on obs time stamp
    ncvar    = tnc.variables
    timename_sim = 'ocean_time'
    utime_sim    = netcdftime.utime(ncvar[timename_sim].units)
    times        = ncvar[timename_sim][:] 
    sdates_sim   = utime_sim.num2date(times)

    #To find closest time step to assimilation time
    sec_assim    = utime_sim.date2num(date_assim)
    diff         = np.abs(sec_assim-times)
    ind          = np.where(diff==diff.min())
    ntime        = pl.array(ind[0]).item()
    print 'ntime > ',ntime,'   date_simul > ',sdates_sim[ind],'   date_assim  > ',date_assim

    z      = tnc.variables['zeta'][ntime,:]

    u_tmp  = tnc.variables ['ubar'][ntime,:]
    u_mask = tnc.variables['mask_u'][:]
    u_tmp[u_mask==0]=0
    #u_tmp[u_tmp.mask]=0

    v_tmp  = tnc.variables['vbar'][ntime,:]
    v_mask = tnc.variables['mask_v'][:]
    v_tmp[v_mask==0]=0    
    #v_tmp[v_tmp.mask]=0

    #calculate u and v in rho points
    u      = pl.zeros_like(z)
    v      = pl.zeros_like(z)

    u[:,1:-1] = (u_tmp[:,:-1]+u_tmp[:,1:])/2.
    u[:,0]    =  u_tmp[:,0]
    u[:,-1]   =  u_tmp[:,-1]

    v[1:-1,:] = (v_tmp[:-1,:]+v_tmp[1:,:])/2.
    v[0   ,:] =  v_tmp [0  ,:]
    v[-1  ,:] =  v_tmp [-1 ,:]


    #to read final mask
    nco       = netCDF4.Dataset(final_grd)
    ncvo      = nco.variables
    mr        = ncvo['mask_rho'][:]
    nco.close()

    #apply the mask
    h = h * mr
    z = z * mr
    u = u * mr
    v = v * mr

    #open directory
    os.system('mkdir -p '+dirout)

    wrt_field(outfile = dirout+'/bathy.bot' ,name =h)  
    wrt_field(outfile = dirout+'/wlevel.bot',name =z)  
    wrt_vec(outfile   = dirout+'/current.bot' ,name1=u,name2=v)  


    ## Write the coordinate file
    # Output File
    outfile=dirout+'/local_rect_coords.bot'
    outf=open(outfile,'w')
    ny,nx=x.shape
    aa=pl.arange(ny-1,-1,-1)
    for iy in aa:
        for ix in range(nx):
            str1=('%.4f' %  x[iy,ix])
            str1=str1+'\n'
            outf.write(str1)

    for iy in aa:
        for ix in range(nx):
            str1=('%.4f' %  y[iy,ix])
            str1=str1+'\n'
            outf.write(str1)  
    outf.close()

###############################################################################################
if real_data:
    # >>>    copy boundary spectra from inner model
    #Find directory of inner to get BOU spectra 
    #inner='/home/shusin3/users/ggarcia/projects/new_river_inlet/swan/201205/16-inner-tide-currents/steady/'
    #inner='/home/shusin3/users/ggarcia/projects/new_river_inlet/swan/201205_ww3_nested/16-inner-tide-currents/steady/'
    inner='/home/shusin3/users/moghimi/coupling/wave_verification/local/01-ggarcia/04-inner-tide-currents-final/03-stationary/'
    wav_base = datetime(2012,05,01,01,0,0);
    units  = 'seconds since '+wav_base.isoformat()[:10]+' '+wav_base.isoformat()[11:]
    utwav  = netcdftime.utime(units)
    timvec = pl.linspace(0,500,501)*3600.0
    wav_dates = utwav.num2date(timvec)
    ind       = pl.array(pl.where((wav_dates>=date_assim))).squeeze()
    ind0      = ind[0]
    
    if '08_forward' in dirin:
        dirout    = dirin.replace('08_forward','08_forward_swan/08_forward')
    else:
        dirout    = dirin.replace('/03_mem_inp/','/04_wav_adj/')
    spcdir    = inner + str(ind0)
    print spcdir
    os.system('cp '+' '+spcdir+'/local20m.spc '+dirout+'/bnd.spc')
    print 'date_assim > ',date_assim,'    taken_date_wave> ',wav_dates[ind[0]],'\n',spcdir

else:
    #global experiment_spc
    #experiment_spc=False
    #spc_file=local_inp+'/obs_wave/syn/wave_limit_test/spc_inner/local20m_hs.5_t14.spc'
    #spc_file=local_inp+'/obs_wave/syn/wave_limit_test/spc_inner/local20m_hs3_t4.spc'
    spc_file=local_inp+'/obs_wave/syn/wave_limit_test/spc_inner/local20m_hs1_t7.spc'
    print '\n\n >>>  experiment_spc is True >>>>> it should be False by default!'
    print '         >>>>>>>> ******* this experimental spc is taking for swan  ****** <<<<<<'
    print '         >>>>>>>> *******  consider  experiment_spc=False  ???      ****** <<<<<<'
    os.system('ln -sf '+spc_file+'  '+dirout+'/bnd.spc' )
    print ' global experiment_spc  \n',spc_file
############################################



#os.system('mv '+' '+dirs[9]+'/tmp.inp '+dirout)
#os.system('mv '+' '+dirs[9]+'/log1.txt '+dirout)

    
