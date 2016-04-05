#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main script for doing assimilation ####
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

#############################################
# Saeed Moghimi; moghimis@gmail.com
# Logs:
# 1.0 03/25/2013 02:14:41 PM    prepared for general use
# 2.0 04/04/2015 09:14:41 AM    complete overhall- adding python assimilation routine
# 2.1 03/09/2016 				[David Honegger] Additional commenting
# 3.0
# 4.0 
#########################################################
# Import python packages to use
import os,sys
import glob
import pylab as pl
import netCDF4
import netcdftime

# Declare global access to some logical variables
global roms_swan
global asim_sar
global asim_wav
global asim_swf
global real_data
global uv_curv
##############################################################
############## Do not change below this line ################# 
######################### Main code ##########################
##############################################################
# Declare global access to more variables
global base_dir
global inp_dir,local_inp
global scr_dir
global mk_bathy    
global adj_bathy 
global do_run
global adj_mem
global do_assim_matlab
global do_pre_assim_sar
global do_pre_assim_wav
global do_assim_python
#######################
# Declare global access to more variables
global mat2prior
global pre_wav_scr 
global adj_wav_scr  
global pre_swf_scr
global ocean_exe
global run_type
global roms_inp
global run_forwards
global logfile
global final_grd
global bry_file

##########################################
# Clean up from previous runs and import from:
# "base_info.py"; "__init__.py"; "assim_module.py"
##########################################
try:
    os.system('rm base_info.pyc'  )
except:
    pass
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
##########################################
try:
    os.system('rm __init__.pyc'  )
except:
    pass

if '__init__' in sys.modules:  
        del(sys.modules["__init__"])
from __init__ import *
##########################################
try:
    os.system('rm assim_module.pyc'  )
except:
    pass

if 'assim_module' in sys.modules:  
        del(sys.modules["assim_module"])
from assim_module import *


##########################################
# Rename imported variables to global counterparts
##########################################
base_dir  = base_info.base_dir
inp_dir   = base_info.inp_dir
scr_dir   = base_info.scr_dir
prior     = base_info.prior
final_grd = base_info.grd
bry_file  = base_info.bry
##########################################
# Rename imported variables to local counterparts
##########################################
###  Control variables
# number of iterations
nstart     = base_info.nstart
nend       = base_info.nend
############  List of available servers  #####
servers    = base_info.servers
###########  Type of operation  (TWIN test or Real data assimilation)   >>>
real_data  = base_info.real_data # False: use synthetic data  True: use read observation
uv_curv    = base_info.uv_curv   # to use curvilinear grid to extract align curve velocity  
if not real_data:
   uv_curv = False      
############  Which model is going to be assimilated >>>>>>>>>>>>>>>>>>>>>>>>
asim_sar   = base_info.asim_sar      # 2. ROMS alone assimilation (waves=False)
asim_wav   = base_info.asim_wav      # 3. SWAN alone assimilation (change the option in
asim_sar2  = base_info.asim_sar2    #In case of assimilating a second SAR image
asim_swf   = base_info.asim_swf     #assimilate swift data 
##############################################
# Make bathymetry parameters global
global equal_space
global N                
global Li              
global Lj          
global Lz           
global sar_err_reduction
global radar_err_reduction
global cov_hh_limit
##############################################
equal_space   = base_info.equal_space
N             = base_info.N     
Li            = base_info.Li    
Lj            = base_info.Lj
Lz            = pl.array(base_info.Lz )
sar_err_reduction   = pl.array(base_info.sar_err_reduction)
radar_err_reduction = pl.array(base_info.radar_err_reduction)
cov_hh_limit  = pl.array(base_info.cov_hh_limit)
########################################################################
### Read in cut window #####
global ijxy

######### Read in script names #############
main_scr    = base_info.main_scr
mk_bathy    = base_info.mk_bathy 
adj_bathy   = base_info.adj_bathy
do_run      = base_info.do_run        
adj_mem     = base_info.adj_mem
############################################
do_pre_assim_sar = base_info.do_pre_assim_sar
do_pre_assim_wav = base_info.do_pre_assim_wav
do_assim_python  = base_info.do_assim_python
###########################################
mat2prior   = base_info.mat2prior
#SWAN
pre_wav_scr = base_info.pre_wav_scr
adj_wav_scr = base_info.adj_wav_scr
#SWIFT
pre_swf_scr = base_info.pre_swf_scr
#
ocean_exe   = base_info.ocean_exe
roms_inp    = base_info.roms_inp
run_type    = base_info.run_type
##########################################
# log file (defined in "__init__.py"); logf() appends to logfile
logf('>>>>>>>>>>>>>> start new assimilation at > ',datetime.datetime.now().isoformat() ,logfile)
logf('base_dir',base_dir ,logfile)
logf('mk_bathy',mk_bathy ,logfile)
logf('adj_bathy',adj_bathy ,logfile)
logf('do_run',do_run ,logfile)
logf('adj_mem',adj_mem ,logfile)
logf('mat2prior',mat2prior ,logfile)
logf('pre_wav_scr',pre_wav_scr ,logfile)
logf('adj_wav_scr',  adj_wav_scr ,logfile)
logf('pre_swf_scr',pre_swf_scr ,logfile)
logf('ocean_exe',ocean_exe ,logfile)
logf('run_type',run_type ,logfile)
logf('roms_inp',roms_inp ,logfile)
#########################################
#copy inp_dir to local_inp
if os.path.exists(local_inp):
    print '>>>>> Local input directory exist ...'
    print '>>>>> Check if OBS files are set right ...'
else:
    print '>>>>> Copying inp_dir to local_inp ...'
    os.system('cp -r ' +inp_dir+'   '+local_inp)

if not real_data:
    uv_curv = False    
####### Main ################
sys.exit()
if __name__ == '__main__':
#main_seq=True
#if main_seq:
    #__main__
    for itr in range(nstart,nend+1):
        logf('**** > Start iteration > ',str(itr)+'  from > '+str(nend) ,logfile)
        # check if we need to create new directory structure
        base=base_dir+'/run_'+str(1000+itr)+'/'
        if os.path.exists(base):
            print 'assimilation directory exist'
            sys.exit('I have to think what to do now !')
        else:
            dirs=create_dirs(itr,base)
        
        # for the first time prior will be copied from inp directory
        if itr==0:
            os.system('cp '+local_inp+'/const/'+prior+' '+dirs[0]+'/prior.nc')
        else:
            #copy new_perior from previous iteration
            new_perior=base_dir+'/run_'+str(1000+itr-1)+'/06_mat2prior/new_prior.nc'
            os.system('cp '+' '+new_perior+' '+dirs[0]+'/prior.nc' )
	
	###### >>>>>>>>  Start ASSIMILATION  <<<<<<<<< #############
	#Asign constants from base_info file
	try:
	    dep_ij=Lz[itr]
	except:
	    dep_ij=Lz.min()
	############################################################
	try:
	    sar_err_reduce   = sar_err_reduction[itr]
	    radar_err_reduce = radar_err_reduction[itr]
	except:
	    sar_err_reduce   = sar_err_reduction.max()
	    radar_err_reduce = radar_err_reduction.max()
	############################################################
	try:
	    cov_hh_limit_ij=cov_hh_limit[itr]
	except:
	    cov_hh_limit_ij=cov_hh_limit.min()            
	############################################################
	logf('>sar_err_reduce   ',str(sar_err_reduce) ,logfile)
	logf('>radar_err_reduce ',str(radar_err_reduce) ,logfile)
	logf('>dep_ij '          ,str(dep_ij) ,logfile)
	logf('>cov_hh_limit_ij ' ,str(cov_hh_limit_ij) ,logfile)
	############################################################
	#Taking care of whether flood or ebb data need to take into account
	#Works only for SYN data now
	############################################################
	if real_data:
	    obs_syn=local_inp+'/obs/sar/'+'uASAR.nc'
	    obs_flood=local_inp+'/obs/sar/'+base_info.obs_flood
	    obs_ebb=local_inp+'/obs/sar/'+base_info.obs_ebb
	else:                
	    obs_syn=local_inp+'/obs/syn/'+'syn1nri_his.nc'
	    obs_flood=local_inp+'/obs/syn/'+base_info.obs_flood
	    obs_ebb=local_inp+'/obs/syn/'+base_info.obs_ebb
	#############################################################    
	try:
	    tide_case=base_info.choose_tide[itr]
	except:
	    tide_case=base_info.choose_tide[-1]

	if  tide_case=='ebb':
	    logf(' ', 'Taking ebb obs file.' ,logfile)
   	    os.system('rm '+obs_syn)
	    os.system('ln -s '+obs_ebb+'  '+obs_syn)
	    logf('assim obs file ',obs_ebb,logfile)
	elif tide_case=='flood':
	    logf(' ', 'Taking flood obs file.' ,logfile)
	    os.system('rm '+obs_syn)
	    os.system('ln -s '+obs_flood+'  '+obs_syn)
	    logf('assim obs file ',obs_flood ,logfile)

	make_bathy(dirs,N,Li,Lj,dep_ij,equal_space)
	#sys.exit()
	bathy_adj(dirs,cov_hh_limit_ij)
	#sys.exit()
	##############################################################
	#do_runs(dirs)
	#do_runs_roms(dirs,servers,nrun=50)
	do_runs_roms(dirs,servers)
	next = check_runs(dirs[3],list=None)
	
	members_adj(dirs,nobs=0)  # nobs is the number of data index in data  nc file (0 for the first time index) 
	do_pre_assim_cur (itr, dirs) 
	
	if asim_sar2:
	    members_adj(dirs,nobs=1)
	
	if asim_swf:
 	    roms2swift(dirs)
	    do_pre_assim_swift(itr, dirs)
	
	if asim_wav:
	    pre_swan_run(dirs)
	    do_swan_run (dirs,servers)
	    #final check over all again (because there might be remaining)
	    next=check_runs(dirs[9],list=None)
	    do_swan_adj(dirs)
	    do_pre_assim_wave(itr, dirs) 

	#Matlab routine 
	#do_assimilate_matlab (dirs,sar_err_reduce,radar_err_reduce,asim_sar,asim_wav,asim_swf)
	#PYTHON routine
	do_assimilate_py(itr, dirs)
	mk_new_prior(dirs)

logf(' > Finish this part ',' With Success < ' ,logfile)

sys.exit() #  
for itr in [0,1,2,3,4,5]:
    print 'Start forward run > ', itr,'  from > ', nend
    base = base_dir+'/run_'+str(1000+itr)+'/'
    dirs = create_dirs(itr,base)
    os.system('cp -f '+' '+dirs[6]+'/new_prior.nc'+' '+dirs[8]+'/grd.nc' )
    do_forward_run(dirs,'tirigan')

#SWAN Forward RUN
for itr in [0,6]:#,4,5,6,7,8]:
    print 'Start SWAN forward run > ', itr,'  from > ', nend
    base = base_dir+'/run_'+str(1000+itr)+'/'
    dirs = create_dirs(itr,base)
    do_forward_swan_run(dirs,'tirigan')

#SWAN Forward RUN adjust and netcdf output
for itr in [0,6]:#,4,5,6,7,8]:
    base = base_dir + '/run_'+str(1000+itr) + '/'
    dirs = create_dirs(itr,base)    
    logf('do_swan_adj()','  > Adjust forward swan run ..' ,logfile)
    sij=str1=' '.join('%.0f' % n for n in ijxy)
    os.system('cd '+dirs[11]+'; echo "'+ sij+ '" &> param.inp  ' )
    os.system('cp -f '+scr_dir+'/'+adj_wav_scr+' '+dirs[11])
    os.system('cp -f '+scr_dir+'/base_info.py    '+dirs[11])
    os.system('cd '+dirs[11]+'; python  '+adj_wav_scr +' '+ dirs[11] +'/08_forward/  >> log_'+adj_wav_scr+'.txt')

logf(' > Finish forward run too ',' With Success < ' ,logfile)
############## END of reg code ###########################################
sys.exit() #  <--- do not remove this
##########################################################################
###  >>>  Here is some useful code snips
### not part of the code
### ready to paste in ipython
clean_runs = True
if clean_runs:
   nstart1=1
   nend1=20 
   for itr in range(nstart1,nend1+1):
        print 'Start cleaning run > ', itr,'  from > ', nend
        base=base_dir+'/run_'+str(1000+itr)+'/'
        dirs=create_dirs(itr,base)
        clean_dirs(dirs)

sys.exit()

#not part of the code
#ready to paste in ipython
itr = 3
real_data   = True      # False: use synthetic data  True: use read observation
asim_sar    = True
asim_wav    = True
asim_swf    = True    #assimilate swift data 
tide_case   = 'ebb'
#tide_case   = 'flood'

############################################################
if real_data:
    obs_syn=local_inp+'/obs/sar/'+'uASAR.nc'
    obs_flood=local_inp+'/obs/sar/'+base_info.obs_flood
    obs_ebb=local_inp+'/obs/sar/'+base_info.obs_ebb
else:                
    obs_syn=local_inp+'/obs/syn/'+'syn1nri_his.nc'
    obs_flood=local_inp+'/obs/syn/'+base_info.obs_flood
    obs_ebb=local_inp+'/obs/syn/'+base_info.obs_ebb
#############################################################    
if  tide_case=='ebb':
    logf(' ', 'Taking ebb obs file.' ,logfile)
    os.system('rm '+obs_syn)
    os.system('ln -s '+obs_ebb+'  '+obs_syn)
    logf('assim obs file ',obs_ebb ,logfile)
elif tide_case=='flood':
    logf(' ', 'Taking flood obs file.' ,logfile)
    os.system('rm '+obs_syn)
    os.system('ln -s '+obs_flood+'  '+obs_syn)
    logf('assim obs file ',obs_flood ,logfile)
#############################################################
base = base_dir+'/run_'+str(1000+itr)+'/'
dirs = create_dirs(itr,base)
#############################################
members_adj(dirs,nobs=0)  
do_pre_assim_cur(itr, dirs)

pre_swan_run(dirs)
do_swan_run (dirs,servers)
next=check_runs(dirs[9],list=None)
do_swan_adj(dirs)
do_pre_assim_wave(itr, dirs)

roms2swift(dirs)
do_pre_assim_swift(itr, dirs)

do_assimilate_py (itr, dirs) 
mk_new_prior(dirs)

##############################################################
#itr = 4
#obs_ebb = local_inp+'/obs/syn/'+base_info.obs_ebb
#os.system('rm '+obs_syn)
#os.system('ln -s '+obs_ebb+'  '+obs_syn)

itr = 0
obs_flood=local_inp+'/obs/syn/'+base_info.obs_flood
os.system('rm '+obs_syn)
os.system('ln -s '+obs_flood+'  '+obs_syn)

base = base_dir+'/run_'+str(1000+itr)+'/'
dirs = create_dirs(itr,base)

do_assimilate_py (itr, dirs) 
mk_new_prior(dirs)

############################################################
