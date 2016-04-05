from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 This script intended to be the main assimilation script
    NEED TO CHECK:
    > 1. Division by numbers 4/5 check not to divide 
    > 2. check cos and sin for degree or radian
    > 3. 
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "moghimis@gmail.com"
#####################################################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 03/25/2013 02:14:41 PM   
#
#
#
import os,sys
############################################
os.system('rm base_info.pyc'  )
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
###########################################
import matplotlib
matplotlib.use('Agg')

import glob
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import netcdftime
from   scipy.optimize import fsolve
from   datetime import datetime
global g
import pickle
#

args = sys.argv
itr = int(args[1])    

sys.path.append("py/") 

#control vars
LOCALIZE         = True
real_data        = base_info.real_data
assim_wav        = base_info.asim_wav
assim_cur        = base_info.asim_sar
assim_swf        = base_info.asim_swf
localize_len     = base_info.Localize_len
radar_err_reduce = base_info.radar_err_reduction[itr]
sar_err_reduce   = base_info.sar_err_reduction[itr]
swift_err_reduce = base_info.swift_err_reduction[itr]
##
#Set input parameters
#base_dir = '/home/shusin5/users/moghimi/assimilation/assim_local/real_data_10days_v2/01_main_ebb/'
#inp_dir  = base_dir + '/inp/'
#scr_dir  = '/home/server/pi/homes/moghimi/work/00-projs/01-muri/00-progs/cowast/project2/07-assim/scr_generations/set4_all/pypysim/'
#
base_dir  = base_info.base_dir
scr_dir   = base_info.scr_dir
prior     = base_info.prior
final_grd = base_info.grd
#                   
######################################################
####################################################
#Set input parameters
#pysim_inp =  base_dir + '/pysim_inp.txt'
#fdata1=open(pysim_inp)
#for  line in fdata1.readlines():
#    if 'itr'    in line: itr          = int  (line.split()[-1]) 
#fdata1.close()
#
inp_dir  = base_dir + '/inp/'
run_id   = '/run_'+str(1000+itr) 
assimilate_dir = base_dir+run_id+'/05_assimilate/'
#####################################################
#import calc_k
g = 9.8126
#sys.path.append(os.environ["INP_DIR"]+'/scr/py')
class component:
    def __init__(self,num):
        self.x=np.zeros((num),dtype='float')
        self.y=np.zeros((num),dtype='float')
        self.data=np.zeros((num),dtype='float')
        self.s=np.zeros((num),dtype='float')
#
# lump observations together
class comp_list:
    def __init__(self):
        self.x=    np.zeros((0)    ,dtype='float')
        self.y=    np.zeros((0)    ,dtype='float')
        self.data= np.zeros((0)    ,dtype='float')
        self.s=    np.zeros((0)    ,dtype='float')
        self.f=    np.zeros((0)    ,dtype='float')
#
# compute model covariances
#
def mycov(a,b):
    """
    Cab = myCov(a,b)
    
    Computes sample covariance 'a' and 'b'.  Inputs should have dimensions
    as follows...
        size(a) = [ M, N ]
        size(b) = [ P, N ]
    That is, each matrix consists of N samples of an Mx1 (or Px1) vector.
    Note M and P need not be equal (unlike the builtin matlab function cov.m).
    Output 'Cab' is the sample covariance, which has dimensions MxP.
    
    a=2  * np.random.rand(3 ,5)
    b=2.2* np.random.rand(4 ,5)
    """
    if(a.shape[1] != b.shape[1]):
        print 'inputs must have same number of samples (2nd dimension)'
    
    [any,anx] = a.shape
    [bny,bnx] = b.shape

    a_mean = np.tile(a.mean(1),a.shape[1]).reshape (a.shape[1],a.shape[0]).T
    b_mean = np.tile(b.mean(1),b.shape[1]).reshape (b.shape[1],b.shape[0]).T
    
    da = a - a_mean
    db = b - b_mean
    N  = anx
    cab = (1.0/(N-1.0)) * np.dot(da,db.T)

    return cab
#
def omega_vec(Cpp,dist,L):
    """
    Cpp = omega_vec(Cpp,dist,L)
    helper function: decorrelation for long length scales (matrix to be
    multiplied onto prior covariance, element by element
    """ 
    if(L==0): return
     
    a   = np.sqrt(10.0/3.0)*L
    b   = dist
    boa = b/a
    c   = 0*Cpp
   
    ind = np.where ((0 <= b) & (b <= a))
    c[ind] =                         \
          -(1.0/4.0)*boa[ind]**5.0   \
          +(1.0/2.0)*boa[ind]**4.0   \
          +(5.0/8.0)*boa[ind]**3.0   \
          -(5.0/3.0)*boa[ind]**2.0   \
          +1.0
    
    ind = np.where ((a < b) &(b <= 2.0 * a)) 
    c[ind] = (1.0 /12.0)*boa[ind]**5.0      \
             -(1.0/2.0 )*boa[ind]**4.0      \
             +(5.0/8.0 )*boa[ind]**3.0      \
             +(5.0/3.0 )*boa[ind]**2.0      \
             - 5.0*boa[ind]                 \
             + 4.0                          \
             -(2.0/3.0)*boa[ind]**(-1.0);

    Cpp=Cpp*c
    return Cpp                  
#
allmeas = {}
if assim_wav:
    print ' > Assim. WAV'
    print '   > Reduce WAV ', radar_err_reduce

    print '    > ','k'
    #read wave dict
    wav_member_dir = base_dir+run_id+'/04_wav_adj/'
    namep      = 'wav_k.nc'
    wpick_name = wav_member_dir + namep
    wavenc     = netCDF4.Dataset(wpick_name)
    wavv       = wavenc.variables
    wavnum     = len(wavv['k_x'][:])
    wmeas   = component(wavnum) 
    wmeas.data = wavv['k_data'][:]
    wmeas.x    = wavv['k_x'][:]
    wmeas.y    = wavv['k_y'][:]
    wmeas.s    = wavv['k_s'][:] * radar_err_reduce
    wmeas.f    = wavv['k_f'][:]
    wmeas.model = wavv['k_model'][:]
    allmeas['k'] = wmeas

if assim_cur:
    cur_member_dir = base_dir+run_id+'/04_mem_adj'
    print ' > Assim. CUR'
    print '   > Reduce CUR ', sar_err_reduce
    
    for field in ['u','v']:
    #for field in ['v']:
        print '    > ',field
        cpick_name = cur_member_dir+'/cur_'+field+'.nc'
        curnc     = netCDF4.Dataset(cpick_name)
        curv       = curnc.variables
        curnum     = len(curv[field+'_x'][:])
        cmeas   = component(curnum) 
        cmeas.data  = curv[field+'_data'][:]
        cmeas.x     = curv[field+'_x'][:]
        cmeas.y     = curv[field+'_y'][:]
        cmeas.s     = curv[field+'_s'][:] * sar_err_reduce
        cmeas.f     = curv[field+'_f'][:]
        cmeas.model = curv[field+'_model'][:]
        allmeas[field] = cmeas


    if base_info.sar_const_err is not None:
        print '    > SAR const err =', base_info.sar_const_err
        allmeas['u'].s =  base_info.sar_const_err * np.ones_like(allmeas['u'].s)
        allmeas['v'].s =  base_info.sar_const_err * np.ones_like(allmeas['v'].s)

    ####################################################################################
    # #SAR err correction when rad data is close
    # #the idea is to increase sar error to decreas its effects when we have wave data close
    if base_info.increase_sar_err_when_waves and assim_cur and assim_wav:
        print '    > Increase SAR err close to Wav data points  dist=',\
        base_info.wav_cur_data_min_dist,'   Coef= ', base_info.cur_data_err_increase_coef
        for isar in range(len(allmeas['u'].x)):
            dist2  = np.sqrt ( (allmeas['k'].x - allmeas['u'].x[isar])**2+\
                           (allmeas['k'].y - allmeas['u'].y[isar])**2  )
            dist_lim = base_info.wav_cur_data_min_dist
            coef     = base_info.cur_data_err_increase_coef
            dist2_min = max(dist2.min(),5)
            if dist2_min < dist_lim :
               allmeas['u'].s[isar] = allmeas['u'].s[isar] * coef * dist_lim  / dist2_min
               allmeas['v'].s[isar] = allmeas['v'].s[isar] * coef * dist_lim  / dist2_min
####################################################################################
###
cur_member_dir = base_dir+run_id+'/04_mem_adj/'
namep          = 'cur_members_prior.nc'
cpick_name     = cur_member_dir + namep
curnc          = netCDF4.Dataset(cpick_name)
curv           = curnc.variables
hc             = curv['h_mems'][:]
xc             = curv['x_rho' ][:]
yc             = curv['y_rho' ][:]
###
if assim_swf:
    swf_member_dir = base_dir+run_id+'/04_swf_adj'
    print ' > Assim. SWF'
    print '   > Reduce SWF ', swift_err_reduce  
    for field in ['us','vs']:
    #for field in ['v']:
        print '    > ',field
        cpick_name = swf_member_dir+'/swf_'+field+'.nc'
        curnc     = netCDF4.Dataset(cpick_name)
        curv       = curnc.variables
        curnum     = len(curv[field+'_x'][:])
        cmeas   = component(curnum) 
        cmeas.data  = curv[field+'_data'][:]
        cmeas.x     = curv[field+'_x'][:]
        cmeas.y     = curv[field+'_y'][:]
        cmeas.s     = curv[field+'_s'][:] * swift_err_reduce
        cmeas.f     = curv[field+'_f'][:]
        cmeas.model = curv[field+'_model'][:]
        allmeas[field] = cmeas

####################################################################################
# #SWIFT err correction when SAR data is close
# #the idea is to increase sar error to decreas its effects when we have wave data close
if base_info.increase_swf_err_when_sar and assim_cur and assim_swf:
    print '    > Increase SWIFT err close to SAR data points  dist=',\
     base_info.swf_cur_data_min_dist,'   Coef= ', base_info.swf_data_err_increase_coef
    for iswf in range(len(allmeas['us'].x)):
        dist2  = np.sqrt ( (allmeas['u'].x - allmeas['us'].x[iswf])**2+\
                           (allmeas['u'].y - allmeas['us'].y[iswf])**2  )
        dist_lim = base_info.swf_cur_data_min_dist
        coef     = base_info.swf_data_err_increase_coef
        
        dist2_min = max(np.abs(dist2.min()),5)
        if dist2_min < dist_lim :
             allmeas['us'].s[iswf] = allmeas['us'].s[iswf] * coef * dist_lim  / dist2_min
             allmeas['vs'].s[iswf] = allmeas['vs'].s[iswf] * coef * dist_lim  / dist2_min
####################################################################################



#nc_prior  = netCDF4.Dataset(inp_dir+'/const/'+prior)
#ncv_prior = nc_prior.variables
#xc     = np.squeeze(ncv_prior['x_rho']    [:])
#yc     = np.squeeze(ncv_prior['y_rho']    [:])
#hc     = np.squeeze(ncv_prior['h']   [:])
#nc_prior.close()
###
###
meas2 = comp_list()
for field in allmeas.keys():
    #print field
    obj=allmeas[field]
    meas2.x    = np.hstack((meas2.x    , obj.x))
    meas2.y    = np.hstack((meas2.y    , obj.y))
    meas2.f    = np.hstack((meas2.f    , obj.f))
    meas2.s    = np.hstack((meas2.s    , obj.s))
    meas2.data = np.hstack((meas2.data , obj.data))
[nn] = meas2.data.shape

#plt.figure()
#plt.scatter(meas2.x,meas2.y,s=50,c=meas2.data,lw=0)   
# Read members depth depth
# Assuming ROMS model area has the biggest coverage
[nyc,nxc,Nc]   = hc.shape

#### construct model for all obs in dicts
meas2_model = np.zeros((nn,Nc)    ,dtype='float')
for il in range(Nc):
    test = np.zeros((0)    ,dtype='float')
    #print il
    for field in allmeas.keys():
        #print field
        obj  = allmeas[field]
        inp  = obj.model[:,il]
        test = np.hstack((test,inp))
    meas2_model[:,il] = test
  
meas2.model = meas2_model
measf       = meas2

#sys.exit()
#-----------------------------------------
# assimilate
#-----------------------------------------   
# We need some measures for observation errors
Cdd    = np.diag(measf.s**2.0)
hvec   = hc.reshape(nxc*nyc, Nc);
model  = measf.model

print ' > Compute model covariances'
Chv    = mycov(hvec ,model)
LCvvL  = mycov(model,model)
# 
xcf    = xc.flatten()
ycf    = yc.flatten()
ng     = len(xcf)
n_meas =len(measf.x)
# 
# 
if LOCALIZE:
    L      = localize_len;
    print ' > Localizing covariances L=', L, 'm'
    if False:
        from omegaf_vec import  omegaf_vec
        Chv  ,distg = omegaf_vec( Chv  , xcf     , ycf     , measf.x, measf.y, L, ng    , n_meas )
        Cdd  ,distx = omegaf_vec( Cdd  , measf.x, measf.y, measf.x, measf.y, L, n_meas, n_meas )
        LCvvL,distx = omegaf_vec( LCvvL, measf.x, measf.y, measf.x, measf.y, L, n_meas, n_meas )
    
    #Chv3   = omega_dist_vec( Chv,  xcf    ,     ycf, measf.x, measf.y, L)
    #Cdd3   = omega_dist_vec( Cdd,  measf.x, measf.y, measf.x, measf.y, L )
    #LCvvL3 = omega_dist_vec( LCvvL,measf.x, measf.y, measf.x, measf.y, L )
    # 
    else:
        from distg import distg as dist
        distg = dist(xcf     ,     ycf, measf.x, measf.y, ng    , n_meas )
        distx = dist(measf.x , measf.y, measf.x, measf.y, n_meas, n_meas )
        
        Chv   = omega_vec  (Chv  ,distg,L);
        Cdd   = omega_vec  (Cdd  ,distx,L);
        LCvvL = omega_vec  (LCvvL,distx,L);

###########################################################################
[ndata]=measf.data.shape
noise=np.zeros((ndata,Nc))
for n in range (Nc):
   noise[:,n] = measf.s * np.random.randn(ndata) 

dh=np.zeros_like(hvec)
# assimilate for posterior ensemble.  Add random noise to
# observations to ensure correct posterior ensemble covariance
print ' > Assimilating for posterior ensemble'

ChvinvC  = np.dot( Chv  , np.linalg.inv( LCvvL  + Cdd  ))

for n in range(Nc):
    dh[:,n]= np.dot(ChvinvC,(measf.data+noise[:,n]-model[:,n]))

hpost   = hvec+dh
outh = hpost.mean(1).reshape(nyc,nxc)

if False:
    print ' > Pickle outputs'
    out_pick={'hpost':hpost,'hpri':hvec,'dh':dh,'Chv':Chv, \
              'LCvvL':LCvvL,'Cdd':Cdd, 'measf':measf,\
              'xc':xc, 'yc':yc}
    
    import cPickle as pickle
    pick_name = 'assimilate_out.p'
    pickle.dump( out_pick, open(pick_name , "wb" ) )

###########################################################
hpost = hpost.reshape(nyc,nxc,Nc)
hprio = hvec.reshape (nyc,nxc,Nc)
maskh = (hc <-5.0)

hpost = np.ma.masked_array(hpost,maskh)
hprio = np.ma.masked_array(hprio,maskh)

hpost_stdv = hpost.std(2)
hprio_stdv = hprio.std(2)

file_sufix  = '_real_data-' + str (base_info.real_data)
file_sufix += '_cur-'       + str (base_info.asim_sar)
file_sufix += '_wav-'       + str (base_info.asim_wav)
file_sufix += '_swf-'       + str (base_info.asim_swf)
file_sufix += '_curJ-'      + str (base_info.jump_cur)
file_sufix += '_wavJ-'      + str (base_info.jump_wav) 
file_sufix += '_CurErr-'    + str (sar_err_reduce) 
file_sufix += '_WavErr-'    + str (radar_err_reduce) 
file_sufix += '_SwfErr-'    + str (swift_err_reduce)

file_sufix += '_L-'         + str (localize_len) 

namep      = 'assimilate_out'+file_sufix+'.nc'
out_name   = assimilate_dir + namep
outnc      = netCDF4.Dataset(out_name,'w',format='NETCDF3_CLASSIC')

outnc.createDimension('nx'    ,nxc   )
outnc.createDimension('ny'    ,nyc   )
outnc.createDimension('nmem'  ,Nc )
outnc.createDimension('nobs'  ,len(measf.data))


p0 = outnc.createVariable('x_rho','f8',('ny','nx',))
p0.missing_value = -9999.0
p0[:]            = xc

p1 = outnc.createVariable('y_rho','f8',('ny','nx',))
p1.missing_value = -9999.0
p1[:]            = yc


p2 = outnc.createVariable('h_post','f8',('ny','nx','nmem'))
p2.missing_value = -9999.0
p2[:]            = hpost

p3 = outnc.createVariable('h_prior','f8',('ny','nx','nmem'))
p3.missing_value = -9999.0
p3[:]            = hprio

p4 = outnc.createVariable('h_post_std','f8',('ny','nx'))
p4.missing_value = -9999.0
p4[:]            = hpost_stdv

p5 = outnc.createVariable('h_prio_std','f8',('ny','nx'))
p5.missing_value = -9999.0
p5[:]            = hprio_stdv

p20 = outnc.createVariable('obs_x','f8',('nobs',))
p20.missing_value = -9999.0
p20[:]            = measf.x

p21 = outnc.createVariable('obs_y','f8',('nobs',))
p21.missing_value = -9999.0
p21[:]            = measf.y


p22 = outnc.createVariable('obs_s','f8',('nobs',))
p22.missing_value = -9999.0
p22[:]            = measf.s


p23 = outnc.createVariable('obs_data','f8',('nobs',))
p23.missing_value = -9999.0
p23[:]            = measf.data


p24 = outnc.createVariable('obs_f','f8',('nobs',))
p24.missing_value = -9999.0
p24[:]            = measf.f

p25 = outnc.createVariable('obs_model','f8',('nobs','nmem',))
p25.missing_value = -9999.0
p25[:]            = measf.model

readme =  ' \n localiz_length='        +  str (localize_len) 
readme += ' \n member_num='            +  str(Nc)
readme += ' \n base_dir     ='         +  base_dir            
readme += ' \n itr     ='              +  str(itr)            
readme += ' \n inp_dir     ='          +  inp_dir            
readme += ' \n real_data='             +  str(base_info.real_data)
readme += ' \n current_opt ='          +  str (base_info.asim_sar)  
readme += ' \n wave_model='            +  str (base_info.asim_wav)
readme += ' \n roms2       ='          +  str (base_info.asim_sar2)
readme += ' \n sar_err_reduction='     +  str (base_info.sar_err_reduction)
readme += ' \n radar_err_reduction='   +  str (base_info.radar_err_reduction) 
readme += ' \n swift='                 +  str (base_info.asim_swf)
readme += ' \n roms1_data_jump='       +  str (base_info.jump_cur)
readme += ' \n wave_data_jump='        +  str (base_info.jump_wav) 
readme += ' \n curve_grid4uv='         +  str (base_info.uv_curv)
readme += ' \n SAR_err_reduce='        +  str (sar_err_reduce) 
readme += ' \n WAV_err_reduce='        +  str (radar_err_reduce) 
readme += ' \n SWF_err_reduce='        +  str (swift_err_reduce)

outnc.history  = 'moghimis@gmail.com  wave data structure for assimilation '+datetime.now().isoformat() + readme
#outnc.att=  infile[-12:]
outnc.close()

args = sys.argv
scr_name = args[0]    
scr_dir1 = os.getcwd()
os.system('cp -fr  ' + scr_name + '    ' + assimilate_dir)
os.system('cp -fr  base_info.py        ' + assimilate_dir)

if True:
    figname = assimilate_dir + '/pic_'+file_sufix
    figname1 = figname+'.png'
    figname2 = figname+'_prior.png'
    ##### posterior
    plt.figure()
    plt.pcolor(xc,yc,outh)
    plt.clim(-1,7)
    plt.colorbar()
    plt.contour(xc,yc,outh,colors='k',levels=np.linspace(-2, 8, 12))
    plt.savefig(figname1,dpi=450)
    ##### prior
    plt.figure()
    plt.pcolor(xc,yc,hvec.mean(1).reshape(nyc,nxc))
    plt.clim(-1,7)
    plt.colorbar()
    plt.contour(xc,yc,hvec.mean(1).reshape(nyc,nxc),colors='k',levels=np.linspace(-2, 8, 12))
    plt.savefig(figname2,dpi=450)

print '  END > '






















# ################################################################3
# if False:
#     xc1 = [0. , 100. , 200., 400.]
#     
#     yc1 = 2 * xc1
#     xc2, yc2 = np.meshgrid(xc1, yc1)
#     xc2f = xc2.flatten()
#     yc2f = yc2.flatten()
#     
#     
#     mod1 = (xc2f + yc2f + 100)/(xc2f+0.5* yc2f+100)
#     
#     mod = tile ( (mod1 + np.random.randn(len(mod1)) * mod1.mean()),5).reshape (len(mod1),5) 
#     
#     LCvvL2 = mycov(mod,mod)
#     
#     plt.figure()
#     plt.pcolor(LCvvL2)
#     plt.colorbar()
#     
#     L=1000
#     n_meas = len(xc2f)
#     LCvvL3,distx = omegaf_vec( LCvvL2, xc2f , yc2f , xc2f , yc2f , L , n_meas, n_meas )
#     plt.figure()
#     plt.pcolor(LCvvL3)
#     plt.colorbar()
# 
# 
# if False:
#     plt.figure()
#     plt.imshow(np.flipud(distg))
#     plt.title('distg')
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(distx))
#     plt.title('distx')
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(Chv[::2,::2]))
#     plt.title('Chv')
#     plt.clim(-0.0001,0.0001)
#     plt.colorbar()
#     
#     
#     
#     plt.figure()
#     plt.imshow(np.flipud(Chv3))
#     plt.title('Chv loc')
#     plt.clim(-0.5,0.5)
#     plt.colorbar()
#     
#     
#     plt.figure()
#     plt.imshow(np.flipud(LCvvL))
#     plt.title('LCvvL')
#     plt.clim(-0.0001,0.0001)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(LCvvL2))
#     plt.title('LCvvL locf')
#     plt.clim(-0.0001,0.0001)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(LCvvL2))
#     plt.title('LCvvL locp')
#     plt.clim(-0.0001,0.0001)
#     plt.colorbar()
#     #close('all')
#     
#     #if False:
#     plt.figure()
#     plt.imshow(np.flipud(ChvinvC))
#     plt.title('ChvinvC')
#     plt.clim(-1,1)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(ChvinvC2))
#     plt.title('ChvinvC loc')
#     plt.clim(-1,1)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(( LCvvL  + Cdd  )))
#     plt.title('( LCvvL  + Cdd  )')
#     plt.clim(-0.01,0.01)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(inv_LCvvL_Cdd))
#     plt.title('inv(LCvvL+Cdd) ')
#     plt.clim(-1000,1000)
#     plt.colorbar()
#     
#     plt.figure()
#     plt.imshow(np.flipud(inv_LCvvL_Cdd2))
#     plt.title('inv(LCvvL+Cdd) loc')
#     plt.clim(-1000,1000)
#     plt.colorbar()
