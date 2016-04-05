#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main functions for handling assimilation
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

import sys
import glob
##########################################
# Clean up old runs and import/run base_info.py
try:
    os.system('rm base_info.pyc'  )
except:
    pass
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
from base_info import *
##########################################import glob
import numpy as np
from __init__ import *

#############################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0  Only modules here 
# 2.0  Try to got read of args as much as possible
# 3.0   
# 4.0  
#########################################################
############ Functions #######################
#print str1
#############################################

def create_dirs(itr,base): # Creates output directories for each assimilation iteration
    logf('create_dirs()',  base ,logfile)
    dirs=[]
    dirs.append(base+ '00_prior')
    dirs.append(base+ '01_bat_inp')
    dirs.append(base+ '02_bat_adj')
    dirs.append(base+ '03_mem_inp')
    dirs.append(base+ '04_mem_adj')
    dirs.append(base+ '05_assimilate')
    dirs.append(base+ '06_mat2prior')
    dirs.append(base+ '07_post')
    dirs.append(base+ '08_forward')
    #wave directories
    dirs.append(base+ '04_wav_adj')
    dirs.append(base+ '04_swf_adj')
    dirs.append(base+ '08_forward_swan')
    
    for idir in(range(len(dirs))):
        os.system('mkdir -p '+ dirs[idir])
    
    # To keep original scripts and base_info in one place
    os.system('mkdir -p ' + base+'/scr')
    os.system('cp   -fr ' +scr_dir+'/*.py  ' +base+'/scr')
    os.system('cp       ' +scr_dir+'/*.m   ' +base+'/scr')
    os.system('cp -fr ' +scr_dir+'/mat     ' +base+'/scr/mat')
    os.system('cp -fr ' +scr_dir+'/py      ' +base+'/scr/py')
    return dirs

def make_bathy(dirs,nmember,len_i,len_j,dep_ij,equal_space): # Run the mk_bathy() script (currently in matlab)
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('make_bathy()','  > run matlab for generating new set of bathy files ..' ,logfile)
    input_file=dirs[0]+'/input_param.txt'
    filein=open(input_file,'w')
    filein.write(scr_dir   +'\n')
    filein.write(local_inp +'\n')  
    filein.write(str(nmember)+'\n')  #'Num_members_N '
    filein.write(str(len_i)  +'\n')  #'Lengthi_Li    '
    filein.write(str(len_j)  +'\n')  #'Lengthj_Lj    '
    filein.write(str(dep_ij) +'\n')  #'Depthz_Lz     '
    if equal_space:  
        filein.write('1\n')     
    else:
        filein.write('0\n')    
    filein.close()
    
    if not equal_space:
        os.system('cp '+local_inp+'/const/nri_curv_grid.nc '+dirs[0]+'/nri_curv_grid.nc')
    os.system('cp ' +scr_dir+'/'+mk_bathy     +' ' +dirs[0])
    #os.system('cp ' +scr_dir+'/mat/startup.m '+' ' +dirs[0])

    if False:
        os.system('cd '+dirs[0]+'; matlab10b -nodesktop -nosplash >logmkbathy1.txt  < '+ mk_bathy)
    else:
        matlab=' matlab '
        comm='ssh irarum  "cd '+ dirs[0]+';'+matlab+' -nodesktop -nosplash >logasim1.txt  < '+mk_bathy +'; touch bathydone "'
        os.system(comm)

def make_pybathy(dirs,nmember,len_i,len_j,dep_ij,equal_space): 
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('make_bathy()','  > run matlab for generating new set of bathy files ..' ,logfile)
    input_file=dirs[0]+'/input_param.txt'
    filein=open(input_file,'w')
    filein.write(scr_dir   +'\n')
    filein.write(local_inp +'\n')  
    filein.write(str(nmember)+'\n')  #'Num_members_N '
    filein.write(str(len_i)  +'\n')  #'Lengthi_Li    '
    filein.write(str(len_j)  +'\n')  #'Lengthj_Lj    '
    filein.write(str(dep_ij) +'\n')  #'Depthz_Lz     '
    if equal_space:  
        filein.write('1\n')     
    else:
        sys.exit('  >>> Pyhon bathy is not accepting curvilinear steretching of bottom features only works for equal_space=True')    
    filein.close()
    
    matlab=' matlab '
    comm='ssh irarum  "cd '+ dirs[0]+';'+matlab+' -nodesktop -nosplash >logasim1.txt  < '+mk_bathy +'; touch bathydone "'
    os.system(comm)








def bathy_adj(dirs,cov_hh_limit_ij):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('bathy_adj()','  > adjust bathy files for ROMS ..' ,logfile)
    os.system('cp '+scr_dir+'/'+adj_bathy+' '+dirs[2])
    inp_dirr = local_inp + '/const'
    tmpgrd   = inp_dirr  + '/'+ base_info.grd                    
  
    sij = str1=' '.join('%.0f' % n for n in ijxy)
    sij = sij +' '+str(cov_hh_limit_ij)
    os.system('cd '+dirs[2]+'; echo "'+ sij+ '" &> param.inp  ' )
    comm='cd '+dirs[2]+'; python  '+adj_bathy+' ' +dirs[1]+' '+dirs[2]+' '+ tmpgrd+ ' &>logbathyadj1.txt'
    print comm
    os.system(comm)

def do_runs_roms(dirs,servers):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_runs_roms()','  > submitting jobs to servers (ROMS runs py scr) ..' ,logfile)
    inp_dirr = local_inp + '/const'
    roms     = inp_dirr  + '/'+ base_info.ocean_exe
    romsin   = inp_dirr  + '/'+ base_info.roms_inp
    brytmp   = inp_dirr  + '/'+ base_info.bry              #dep_min=-1.0
    varinfo  = inp_dirr  + '/'+ base_info.varinfo

    nmem = nall = N-1
   
    if run_type == '3D':
        nrun  = 21
    else:
        nrun  = (N//len(servers))+1

    imem = 0
    #imem = 105
    while(nall>0):
        print nall 
        dir2check=[]
        for server in servers:
            for irun in range(nrun):
                if imem > nmem: return
                msg='    > Start ROMS run for ... '+server+' '+str(irun+1)+' from '+ str(nrun)+ ' total submitted '+ str(imem+1)
                logf('',msg ,logfile)
                
                dirmem=dirs[3]+'/member'+str(1000+imem+1)
                os.system('mkdir -p '+ dirmem)
                #os.system('ln -s '+romsin+'  '+ dirmem+'/param.in' )
                os.system('cp -f '+romsin +'  '+ dirmem+'/param.in' )
                os.system('ln -sf '+brytmp+   '  '+ dirmem+'/bry.nc' )
                os.system('ln -sf '+varinfo+' '+ dirmem+'/varinfo.dat')
                os.system('ln -sf '+ inp_dirr+'/'+ base_info.initial+   '  '+ dirmem+'/initial.nc')

                #link bathy file
                grd=dirs[2]+'/bath'+str(1000+imem+1)+'.nc'
                os.system('ln -s '+grd+   '  '+ dirmem+'/grd.nc' )
                if run_type == '3D':
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.rain+      '  '+ dirmem+'/rain.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.cloud+     '  '+ dirmem+'/cloud.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.Qair+      '  '+ dirmem+'/Qair.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.Uwind+     '  '+ dirmem+'/Uwind.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.Vwind+     '  '+ dirmem+'/Vwind.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.lwrad_down+'  '+ dirmem+'/lwrad_down.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.net_swrad+ '  '+ dirmem+'/net_swrad.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.Tair+      '  '+ dirmem+'/Tair.nc')
                    os.system('ln -sf '+ inp_dirr+'/'+ base_info.Pair+      '  '+ dirmem+'/Pair.nc')
                    os.system('echo '+  server + '&>'+ dirmem+ '/machinefile')
                    #comm = 'ssh '+server+'  "cd '+ dirmem+'; mpirun -np 1 '+roms+'  param.in > runlog ; echo "run comp" &> flagdone "'
                    comm = 'ssh '+server+'  "cd '+ dirmem+'; '+  \
                        ' source  /home/server/local/apps/pgicdk122-centos5/linux86-64/current/mpi.sh; '+ \
                        ' mpirun  -machinefile machinefile -np 1 ' + roms+\
                        ' param.in > runlog ; echo "run comp" &> flagdone "'
                    
                    os.system(comm+'&')   
                    #print comm
                else:
                    comm='ssh '+server+'  "cd '+ dirmem+'; '+roms+' <  param.in > runlog ; echo run_comp > flagdone "'
                    #print comm 
                    os.system(comm+'&')   

                os.system('sleep 1')
                dir2check.append(dirmem)
                imem+=1 
                nall=nall-1
        next=check_runs(dirs[3],list=dir2check)

def do_swan_run(dirs,servers):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_swan_run()','  > Run SWAN members ..' ,logfile)
    dirlist = glob.glob(dirs[9]+'/mem*')
    dirlist.sort()
    inp_dirr = local_inp + '/const'
    swan          = inp_dirr  + '/'+ base_info.wave_exe
    swan_inp_file = inp_dirr  + '/'+ base_info.wave_inp_file
    #    
    #'/home/shusin3/users/moghimi/projects/nri_swn/swan/swan4081/swan.exe'
    nmem=len(dirlist)
    imem=0
    nall=len(dirlist)
    while(nall>0):
        dir2check=[]
        for server in servers:
            for irun in range(25):
                if imem>nmem-1: return
                msg='    > Start SWAN run for ... '+server+' '+str(irun+1)+' from '+\
                                             ' 25 '+ ' total submitted '+ str(imem+1)
                logf('',msg ,logfile)
                os.system('ln -sf '+swan         +' '+ dirlist[imem]+'/swan' )
                os.system('cp -f  '+swan_inp_file+' '+ dirlist[imem]+'/INPUT')
                #
                #if not real_data:
                #    os.system('rm     '              + dirlist[imem]+'/bnd.spc' )
                #    os.system('cp -rf '+spc_file+'  '+ dirlist[imem]+'/bnd.spc' )
                #
                comm='ssh '+server+'  "cd ' + dirlist[imem]+'; ./swan > runlog ; touch flagdone "'
                os.system(comm+'&')    
                #print comm
                os.system('sleep 1')
                dir2check.append(dirlist[imem])
                imem+=1
                nall=nall-1
        next=check_runs(dirs[9],list=dir2check)

def check_runs(dir2check,list):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('check_runs()','  > Wait for runs to finish ..' ,logfile)
    # check if all already started simulations finished or not.
    # it will not take care of not submitted jobs or it will 
    # check for completion of the members which has runlog already.
    memdir=dir2check+'/member1*'
    if list==None:
        dirlist=glob.glob(memdir)
    else:
        dirlist=list
        logf('    > A specified list will be watched .. ',str(len(dirlist)) ,logfile)
    ################################################
    dirlist.sort()
    nm=len(dirlist)
    next=False
    wait_time=30
    total_wait=0
    while(not next):
        next=True
        for im in (range(nm)):
            romslog=dirlist[im]+'/runlog'
            flag=dirlist[im]+'/flagdone'
            if  os.path.isfile(romslog):
                #flag_exist=os.path.isfile(flag)
                flag_exist=os.path.exists(flag)
                next=(next and flag_exist)
        if not next:
            time.sleep(wait_time) # delays for 10 seconds
            total_wait+=wait_time
            if(np.mod(total_wait,600)==0): 
                #dir_walk=os.walk(dir2check)
                os.system('ls -lRrt  '+dir2check+'   &>ls.log')
                logf('       > The total_wait is (Sec)>>',str(total_wait) ,logfile)
    return next

def members_adj(dirs,nobs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('members_adj()','  > adjust members of ROMS ..' ,logfile)
    #nobs is index of field of the observation cdf file will be fetched. 
    sij=str1=' '.join('%.0f' % n for n in ijxy)
    os.system('cd '+dirs[4]+'; echo "'+ sij+ '" &> param.inp  ' )
    
    if uv_curv and real_data:
        uv_flag='curv'
    else:
        uv_flag='no-curv'    
    
    os.system('cp -f ' +scr_dir+'/'+adj_mem+'  '+dirs[4])
    os.system('cp -f ' +scr_dir+'/base_info.py '+dirs[4])

    comm='cd '+dirs[4]+'; python '+adj_mem+' '+dirs[3] +' '+\
               str(nobs)+' '+uv_flag+'  >>logadj1.txt'
    #print comm
    os.system(comm)

def roms2swift(dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('roms2swift()','  > adjust members for SWIFT ...' ,logfile)
    os.system('cp -f ' +scr_dir+'/'+pre_swf_scr+'  '+dirs[10])
    os.system('cp -f ' +scr_dir+'/base_info.py     '+dirs[10])

    os.system('cd '+dirs[10]+'; python  '+pre_swf_scr+'  '+dirs[3] +' >>logsw1.txt')

def do_pre_assim_cur(itr, dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_pre_assim_cur(itr, dirs)','  > CUR is getting ready for assimilation ..' ,logfile)
    os.system('cp -f ' +scr_dir+'/'+do_pre_assim_sar+'  '+dirs[5])
    os.system('cp -f ' +scr_dir+'/base_info.py          '+dirs[5])
    os.system('cd    '+dirs[5]+'; python -u '+ do_pre_assim_sar+'  '+str(itr)+' >>log_pre_assim_cur.txt')

def do_pre_assim_wave(itr, dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_pre_assim_wave(itr, dirs)','  > WAV is getting ready for assimilation ..' ,logfile)
    os.system('cp -f  ' + scr_dir + '/'+do_pre_assim_wav+'   '+dirs[5])
    os.system('cp -f  ' + scr_dir + '/base_info.py           '+dirs[5])
    os.system('cp -rf ' + scr_dir + '/py                     '+dirs[5])
    os.system('cd   '+dirs[5]+'; python -u '+ do_pre_assim_wav+'  '+str(itr)+' >>log_pre_assim_wav.txt')

def do_pre_assim_swift(itr, dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_pre_assim_swf(itr, dirs)','  > WAV is getting ready for assimilation ..' ,logfile)
    os.system('cp -f  ' + scr_dir + '/'+do_pre_assim_swf+'   '+dirs[10])
    os.system('cp -f  ' + scr_dir + '/base_info.py           '+dirs[10])
    os.system('cp -rf ' + scr_dir + '/py                     '+dirs[10])
    os.system('cd      '+dirs[10]+'; python -u '+ do_pre_assim_swf+
              '  '+str(itr)+' >>log_pre_assim_swf.txt')

def do_assimilate_py(itr, dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_assimilate_py()','  > run python assimilation routine ..' ,logfile)
    os.system('cp -rf '+scr_dir+'/py                     ' + dirs[5])
    os.system('cp -f ' +scr_dir+'/'+ do_assim_python +'  ' + dirs[5])
    os.system('cp -f ' +scr_dir+'/base_info.py           ' + dirs[5])
    #####################################################################################################
    #os.system('cd    '+dirs[5]+'; python  '+ do_assim_python+'  '+str(itr)+' >>log_pre_assim_cur.txt')
    comm = 'ssh tirigan  "cd  '+ dirs[5]+'; /home/server/pi/homes/ggarcia/local/enthought/canopy/User/bin/python  -u '+ do_assim_python+'  '+str(itr)+' ; touch assimdone "'
    print comm
    os.system(comm)
    flag=dirs[5]+'/assimdone'
    while(not os.path.isfile(flag)):
        time.sleep(10) # delays for 10 seconds
    return 

def mk_new_prior(dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('mk_new_prior()','  > making new prior getting ready for next iteration  ..' ,logfile)
    sij=str1=' '.join('%.0f' % n for n in ijxy)
    final_grd = base_info.grd
    os.system('cd ' + dirs[6]  +'; echo "'+ sij+ '" &> param.inp  ' )
    os.system('cp ' + local_inp+'/const/'+final_grd+' '  +dirs[6]+'/final_grd.nc')
    os.system('mv ' + dirs[5]  +'/*  '+dirs[6])
    os.system('cp ' + scr_dir  +'/'+mat2prior+'  '+ dirs[6])
    os.system('cd ' + dirs[6]  +'; python  '+mat2prior+' netcdf  &>log_'+mat2prior+'.txt')

def do_forward_run(dirs,server):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    msg='  > Start forward run for ... \n'+ str(dirs[0])
    #sys.exit(msg)
    logf('do_forward_run()',msg ,logfile)
    inp_dirr = local_inp + '/const'
    roms     = inp_dirr  + '/'+ base_info.ocean_exe
    romsin   = inp_dirr  + '/'+ base_info.roms_inp
    brytmp   = inp_dirr  + '/'+ base_info.bry              #dep_min=-1.0
    varinfo  = inp_dirr  + '/'+ base_info.varinfo
    #    
    os.system('cp -rf '+ romsin + ' '+ dirs[8]+'/param.in' )
    os.system('ln -sf '+ brytmp + ' '+ dirs[8]+'/bry.nc' )
    os.system('ln -sf '+ varinfo +' '+ dirs[8]+'/varinfo.dat')
    if run_type == '3D':
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.initial+   '  '+ dirs[8]+'/initial.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.rain+      '  '+ dirs[8]+'/rain.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.cloud+     '  '+ dirs[8]+'/cloud.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.Qair+      '  '+ dirs[8]+'/Qair.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.Uwind+     '  '+ dirs[8]+'/Uwind.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.Vwind+     '  '+ dirs[8]+'/Vwind.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.lwrad_down+'  '+ dirs[8]+'/lwrad_down.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.net_swrad+ '  '+ dirs[8]+'/net_swrad.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.Tair+      '  '+ dirs[8]+'/Tair.nc')
        os.system('ln -sf '+ inp_dirr+'/'+ base_info.Pair+      '  '+ dirs[8]+'/Pair.nc')
        os.system('echo '+  server + '&>'+ dirs[8]+ '/machinefile')
        #comm = 'ssh '+server+'  "cd '+ dirs[8]+'; mpirun -np 1 '+roms+'  param.in > runlog ; echo "run comp" &> flagdone "'
        comm = 'ssh '+server+'  "cd '+ dirs[8]+'; '+  \
            ' source  /home/server/local/apps/pgicdk122-centos5/linux86-64/current/mpi.sh; '+ \
            ' mpirun  -machinefile machinefile -np 1 ' + roms+\
            ' param.in > runlog ; echo "run comp" &> flagdone "'
    else:
        comm='ssh '+server+'  "cd '+ dirs[8]+'; '+roms+' <  param.in > runlog ; echo "run comp" &> flagdone "'
    #comm='ssh '+server+'  "cd '+ dirs[8]+'; rm  runlog  flagdone; '+roms+' <  param.in > runlog ; echo "run comp" > flagdone "'    
    os.system(comm+'&')     

def pre_swan_run(dirs):
    #prepare and run roms2swan
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('pre_swan_run()','  > SWAN pre-processing routine ..' ,logfile)
    os.system('cp -f '+scr_dir+'/'+pre_wav_scr+'    '+dirs[9]+'/'+pre_wav_scr)
    os.system('cp -f '+scr_dir+'/base_info.py  '+dirs[9])
    
    dirlist = glob.glob(dirs[3]+'/mem*')
    dirlist.sort()
    for dir1 in dirlist:
        dir_temp = dir1.replace('/03_mem_inp/','/04_wav_adj/')
        comm = 'mkdir -p  '+dir_temp+' ;  cd '+dirs[9]+'; python -u '+pre_wav_scr+' '+ dir1 +'  &> log_'+pre_wav_scr+'.txt'
        os.system(comm)
        #print comm

def do_swan_adj(dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_swan_adj()','  > Adjust swan members ..' ,logfile)
    
    sij=str1=' '.join('%.0f' % n for n in ijxy)
    os.system('cd '+dirs[9]+'; echo "'+ sij+ '" &> param.inp  ' )
    os.system('cp -f '+scr_dir+'/'+adj_wav_scr+' '+dirs[9])
    os.system('cp -f '+scr_dir+'/base_info.py    '+dirs[9])
    
    dirlist = glob.glob(dirs[9]+'/mem*')
    dirlist.sort()
    for dir1 in dirlist:
        os.system('cd '+dirs[9]+'; python  '+adj_wav_scr +' '+ dir1 +' >> log_'+adj_wav_scr+'.txt')

def do_forward_swan_run(dirs,server):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    msg='  > Start forward swan run for ... \n'+ str(dirs[0])
    #sys.exit(msg)

    os.system('cp -f '+scr_dir+'/'+pre_wav_scr+'    '+dirs[11]+'/'+pre_wav_scr)
    os.system('cp -f '+scr_dir+'/base_info.py       '+dirs[11])
    comm = 'cd '+dirs[11]+'; python -u '+pre_wav_scr+' '+ dirs[8] +'  &> log_'+pre_wav_scr+'.txt'
    os.system(comm)
    
    logf('do_forward_run()',msg ,logfile)
    inp_dirr = local_inp + '/const'
    swan          = inp_dirr  + '/'+ base_info.wave_exe
    swan_inp_file = inp_dirr  + '/'+ base_info.wave_inp_file

    msg='    > Start SWAN forward run for ... ' + dirs[11]
    logf('',msg ,logfile)
    os.system('ln -sf '+swan         +' '+ dirs[11]+'/08_forward/swan' )
    os.system('cp -f  '+swan_inp_file+' '+ dirs[11]+'/08_forward/INPUT')
    comm='ssh '+server+'  "cd ' + dirs[11]+'/08_forward/; ./swan > runlog ; touch flagdone "'
    os.system(comm+'&')    
    
    #logf('do_swan_adj()','  > Adjust forward swan run ..' ,logfile)
    #sij=str1=' '.join('%.0f' % n for n in ijxy)
    #os.system('cd '+dirs[11]+'; echo "'+ sij+ '" &> param.inp  ' )
    #os.system('cp -f '+scr_dir+'/'+adj_wav_scr+' '+dirs[11])
    #os.system('cp -f '+scr_dir+'/base_info.py    '+dirs[11])
    #os.system('cd '+dirs[11]+'; python  '+adj_wav_scr +' '+ dirs[11] +'/08_forward/  >> log_'+adj_wav_scr+'.txt')

def clean_dirs(dirs):
    logf('>>>>> ',datetime.datetime.now().isoformat() ,logfile)
    logf('do_clean()','  > Cleaning members ..' ,logfile)    
    #Del adjusted bathys
    os.system('cd '+dirs[2]+';  rm bat*')
    #Del roms outputs
    dirlist=glob.glob(dirs[3]+'/mem*')
    for dir1 in dirlist:
        os.system('cd '+dir1+'; rm nri*')
    
    #Del and zip swan outputs
    dirlist=glob.glob(dirs[9]+'/mem*')
    for dir1 in dirlist:
        os.system('cd '+dir1+'; rm bnd.spc')
        os.system('cd '+dir1+'; rm swan')
        os.system('cd '+dir1+'; gzip *')
