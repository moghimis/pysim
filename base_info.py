#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration file for assimilation
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2014, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

### 2D   >>>>>
base_dir = '/home/shusin3/shared/assim/dah_test2_saeed2'
run_type    = '2D'

################################################################################
inp_dir='/home/nmg/nmg/assimilation/set6_python_test_case/inp'
scr_dir='/home/nmg/nmg/assimilation/set6_python_test_case/pysimHonegger'
################################################################################
###################### BASIC DATA  #############################################
# number of iterations
nstart = 0 
nend   = 1
############  List of available servers  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
servers=[#'irarum',
         #'rimus', 
         #'nanum', 
         'tirigan', 
         #'criss',
         #'stanley',
         #'dudu', 
         #'gudea', 
         #'hablum', 
         #'kurum',
         ]

###########  Type of operation  (TWIN test or Real data assimilation)   >>>
real_data  = False    # False: use synthetic data  True: use read observation
uv_curv    = False     # always false <<<<<<<<  to use curvilinear grid to exteract align curve velocity  
############  Which model is going to be assimilated >>>>>>>>>>>>>>>>>>>>>>>>
#Only one of the below options could be TRUE
#Here starts different arrangements
asim_sar   = True     # 2. ROMS alone assimilation (waves=False)
asim_wav   = True     # 3. SWAN alone assimilation (change the option in
                      #    matlab code to not include velocity in assimilation
                      #    but ROMS still needs to run for getting current into 
                      #    SWAN
############  In case of assimilating a SWIFT data        >>>>>>>>>>>>>>>>>>>>
asim_swf  =  False    #assimilate swift data 
############  In case of assimilating a second SAR image  >>>>>>>>>>>>>>>>>>>>
asim_sar2 = False    # always False <<<<<< not working at this point
################################################################################

######  Name of the scripts will be used by main script to treat the files  ####
################################################################################
#scripts name
#ROMS
main_scr          = '00-run_assim_noIter.py'
mk_bathy          = '01-make_bathy_bumps.m'  
adj_bathy         = '02-adjust_bathy_dh.py'
do_run            = '03-doRuns.sh'           
adj_mem           = '04-simplify_roms.py'
########################################
do_assim_python   = '05-py-assimilate.py' 
do_pre_assim_sar  = '05-pre_asim_cur.py'
do_pre_assim_wav  = '05-pre_asim_wav.py'
do_pre_assim_swf  = '05-pre_asim_swf.py'
#######################################
mat2prior   ='06-posterior2prior.py'
#SWAN
pre_wav_scr ='04-roms2swan.py'
adj_wav_scr ='04-simplify_swan.py'
#SWIFT
pre_swf_scr ='04-roms2swift.py'

################################################################################
####  Generating ensemble members (bathymetry)   >>>>>>>>>>>>>>>>>>>>
#making members bathymetry
N = 10     #number of members
###########################################################################
equal_space = True
# perturbation scales
if equal_space:
    #Radius equal to LixLj times number of grid point perturbations    
    Li = 300;      
    Lj = 300;
else:
    Li = 25;      
    Lj = 25;

# depth of bumps
Lz            = 0.5,0.3,0.2,
cov_hh_limit  = 1.5,1.5,1.5,                                # consistent with ebb case

jump_wav = 2  #jump over data points (every other data points for j=2)
jump_cur = 2 
Localize_len = 400   # in meter

###########
increase_sar_err_when_waves = False
if increase_sar_err_when_waves:
    wav_cur_data_min_dist = 40  #[m]
    cur_data_err_increase_coef = 5
########################################
increase_swf_err_when_sar = False
if increase_swf_err_when_sar:
    swf_cur_data_min_dist = 5  #[m]
    swf_data_err_increase_coef = 5
########################################
sar_const_err    = 0.2   #[ms-1]
########################################
## Define assimilation window ###
# [0,-1] is full domain (default)
ix=0,-1  
jy=0,-1

if real_data:
    ########## Real STUFF  ######################################
    # which data set to assimilate in each assimilation iteration
    obs_flood = 'XXXX20120513_132446_stitched_ks.nc'
    obs_ebb   = '20120510_182231_stitched.nc'
    ##
    sar_err_reduction   = 1.0,1.0,1.0,1.0
    radar_err_reduction = 1.0,1.0,1.0,1.0
    swift_err_reduction = 0.05,0.05,0.05,0.05
else:
    ########## SYN STUFF  ######################################
    obs_flood = 'obs_syn1nri_flood.nc'   
    obs_ebb   = 'obs_syn1nri_ebb.nc'
    sar_err_reduction   = 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
    radar_err_reduction = 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
    swift_err_reduction = 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05   # this makes swift error 5 cm/s

    ############################################################
    

e='ebb'
f='flood'
#choose_tide = f,f,f,f,f,f
#choose_tide  = e,f,e,f,e,f
choose_tide  = e,e,e,e,e,e

## Wave model info
wave_exe      = '/swan/' + 'swan_ser'
wave_inp_file = '/swan/' + 'input_rect_local.swn'

## ROMS info
# Forcing files
prior = 'new_prior.nc'
if run_type == '3D':
    ocean_exe = '/roms_3d/' + '30_roms'
    roms_inp  = '/roms_3d/' + 'romsInput_3d.in'
    grd       = '/roms_3d/' + 'nri_regional_grd_roms_20120510_not_smooth.nc'
    bry       = '/roms_3d/' + 'bry_local_v0_3d_from_wl_obs_ver0.nc'
    initial   = '/roms_3d/' + 'clm_init_nri_his_0009_last_ini.nc'
    varinfo   = '/roms_3d/' + 'varinfo.dat'
    ##### METEO stuff ##
    meteo_dir = '/roms_3d/' + '/meteo/'
    cloud = meteo_dir + 'cloud.nc'
    Qair  = meteo_dir + 'Qair.nc'
    Uwind = meteo_dir + 'Uwind.nc'
    Vwind = meteo_dir + 'Vwind.nc'
    lwrad_down = meteo_dir + 'lwrad_down.nc'
    net_swrad  = meteo_dir + 'net_swrad0.8.nc'
    rain  = meteo_dir + 'rain.nc'
    Tair  = meteo_dir + 'Tair.nc'
    Pair  = meteo_dir + 'Pair.nc'
    #
    hisfile_name =  'nri_his_0001.nc'    
    ################
else:
    ocean_exe = '/roms_2d/' + 'coawstS'
    roms_inp  = '/roms_2d/' + 'romsInput_3d_4_2d.in'
    grd       = '/roms_2d/' + 'grd_local_20120510.nc'
    bry       = '/roms_2d/' + 'bry_local_v0_3d_from_wl_obs_ver0.nc'
    initial   = '/roms_2d/' + 'clm_init_nri_his_0009_last_ini.nc'
    varinfo   = '/roms_2d/' + 'varinfo.dat'
    #
    hisfile_name =  'nri_his_0001.nc'    
