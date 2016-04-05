from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pre-process Radar data and put them in vector form.
 - A little bit of cleaning as well
 - Finally interpolated Radar and SWAN will be save in a netcdf file to read by assimilation routine 
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
#####################################################################
import os,sys
import glob
import numpy as np
import netCDF4
import multiprocessing
from   datetime import datetime
import netcdftime
from   scipy.optimize import fsolve

args = sys.argv
itr = int(args[1])    

sys.path.append("py/") 
from   dispersion import *
#
os.system('rm base_info.pyc'  )
if 'base_info' in sys.modules:  
    del(sys.modules["base_info"])
import base_info
####################################################
base_dir  = base_info.base_dir
inp_dir   = base_info.inp_dir
scr_dir   = base_info.scr_dir
prior     = base_info.prior
final_grd = base_info.grd
real_wave = base_info.real_data
jumpp     = base_info.jump_wav
####################################################
#pysim_inp =  base_dir + '/pysim_inp.txt'
#Set input parameters
#fdata1 = open(pysim_inp)
#for  line in fdata1.readlines():
#    print line
#    if 'itr'    in line: itr          = int(line.split()[-1]) 
#fdata1.close()
#####################################################
#import calc_k
g        = 9.8126
run_id   = '/run_'+str(1000+itr)
inp_dir  = base_dir + '/inp/'
#### Funcs
##   interpolation
methodi='csa'
if methodi=='csa':
    import octant.csa as csa
elif methodi=='tri':
    from delaunay import triangulate
#
def interpg(x_old,y_old,data_old,x_new,y_new):
    if methodi=='csa':
        csa_interp = csa.CSA(x_old,y_old,data_old)
        data_new = csa_interp(x_new,y_new)
        return data_new
    elif methodi=='grd':
        from matplotlib.mlab import griddata
        data_new = griddata(x_old.flatten(),y_old.flatten(),data_old.flatten(),x_new,y_new)
    elif methodi=='tri':
        tri=triangulate.Triangulation(x_old,y_old)
        interp_b=tri.nn_extrapolator(data_old)
        data_new = interp_b(x_new,y_new)
        return data_new
#
class component:
    def __init__(self,num):
        self.x=np.zeros((num),dtype='float')
        self.y=np.zeros((num),dtype='float')
        self.data=np.zeros((num),dtype='float')
        self.s=np.zeros((num),dtype='float')
#
#
#### End of Funcs
wav_member_dir = base_dir+run_id+'/04_wav_adj'
data_dir       = inp_dir+'/'

# #### Pickle name and check if exist !
namep      = '/wave_jump-'  + str(jumpp)+'__itr-'+str(itr)+'.p'
pick_name  = wav_member_dir + namep
# if  os.path.isfile(pick_name):
#     sys.exit('Wave file exist !')
###################################
print 'Loading SWAN members '
###################################
memdir  = wav_member_dir+'/member1*'
dirlist = np.array(glob.glob(memdir))
dirlist.sort()
nmemw   = len(dirlist)
#create variables to be read in
# read SIZE from a sample
ncf   = dirlist[0]+'/swout.nc'
nc    = netCDF4.Dataset(ncf)
ncv   = nc.variables
h     = ncv['dep'][:]
xw    = ncv['x'][:];xw1=xw[1,:]
yw    = ncv['y'][:];yw1=yw[:,1]
[ny,nx]=np.squeeze(h).shape
nc.close()
#
uw    = np.zeros((ny,nx,nmemw),dtype='float')
vw    = np.zeros_like(uw)
dw    = np.zeros_like(uw)
Hw    = np.zeros_like(uw)
mw    = np.zeros_like(uw)
dirw  = np.zeros_like(uw)
#
failind=[]
for i in range(nmemw):
    if np.mod(i,50)==0: print ' > WAV member read in ', (i*1.0/nmemw)*100.0, '%'
    try:
        ncf       = dirlist[i]+'/swout.nc'
        nc        = netCDF4.Dataset(ncf)
        ncv       = nc.variables
        uw[:,:,i] = np.squeeze(ncv['velx'][:,:])
        vw[:,:,i] = np.squeeze(ncv['vely'][:,:])
        dw[:,:,i] = np.squeeze(ncv['dep'][:,:])
        mw[:,:,i] = np.squeeze(ncv['mask'][:,:])
        Hw[:,:,i] = np.squeeze(ncv['hsig'][:,:])
        dirw[:,:,i]=np.squeeze(ncv['dir'][:,:])
        nc.close()
    except:
        failind.append(i)
# In case some members failed for whatever reason: discard them
if(len(failind) > 0):
    print 'WARNING WARNING WARNING: the following members failed:'
    for i in failind:
        print 'member -> ', dirlist[i]
    #
    keepind=setdiff1d( range(nmemw),failind)
    uw=uw[:,:,keepind]
    vw=vw[:,:,keepind]
    mw=mw[:,:,keepind]
    dw=dw[:,:,keepind]
    Hw=Hw[:,:,keepind]
    dirw=dirw[:,:,keepind]
    dirlistw=dirlist[keepind]
#mask out land points
uw = np.ma.masked_where(mw==0,uw)
vw = np.ma.masked_where(mw==0,vw)
dw = np.ma.masked_where(mw==0,dw)
Hw = np.ma.masked_where(mw==0,Hw)
dirw=np.ma.masked_where(mw==0,dirw)
maskm = mw[:,:,0]==0
[nyw,nxw,Nw]=dw.shape

print '  > Number of wave members >', Nw
#sys.exit('asdasdasdasda')
######################################
# Read Radar  data
# load the data grid for currents 

if real_wave:
    ###Real OBS
    ncf_wdata = data_dir+'obs_wave/rad/cBathy.nc'
    print 'real data obs file name > ', ncf_wdata
    nc_wdata=netCDF4.Dataset(ncf_wdata)
    ncv_wdata=nc_wdata.variables
    xwd    = np.squeeze(ncv_wdata['x']    [:])
    ywd    = np.squeeze(ncv_wdata['y']    [:])
    fdwd   = np.squeeze(ncv_wdata['fb']   [:])
    kdwd   = np.squeeze(ncv_wdata['k']    [:])
    dirwd  = np.squeeze(ncv_wdata['dir']  [:])
    kd_err = np.squeeze(ncv_wdata['k_err'][:])
    nc_wdata.close()
    ###    
    if jumpp > 1:
        kxx     = jumpp
        kyy     = jumpp
        xwd     = xwd   [  ::kyy,::kxx]
        ywd     = ywd   [  ::kyy,::kxx]
        fdwd    = fdwd  [:,::kyy,::kxx]
        kdwd    = kdwd  [:,::kyy,::kxx]
        dirwd   = dirwd [:,::kyy,::kxx]
        kd_err  = kd_err[:,::kyy,::kxx]   
    ###    
    nkw,nyw,nxw = fdwd.shape
    try:
        mtmp    = fdwd.mask
    except:
        fdwd    = np.ma.masked_where(fdwd > 10, fdwd)
    ###
    mmax    = np.zeros((nyw,nxw))
    kdmax   = np.zeros((nyw,nxw))
    fdmax   = np.zeros((nyw,nxw))
    dirmax  = np.zeros((nyw,nxw))
    kerrmax = np.zeros((nyw,nxw))
    ###
    for ii in range(nxw):
        for jj in range(nyw):
            if not (fdwd[:,jj,ii].mask.all==True):
               indw           = np.argmax(fdwd[:,jj,ii])
               fdmax  [jj,ii] = fdwd [indw,jj,ii]
               kdmax  [jj,ii] = kdwd [indw,jj,ii]
               dirmax [jj,ii] = dirwd[indw,jj,ii]
               kerrmax[jj,ii] = kd_err[indw,jj,ii]
    ###
    #Eliminate data points out side of the model area
    fdmax[ (xwd<xw.min()) | (xwd > xw.max()) | (ywd<yw.min()) |(ywd > yw.max()) ]=np.NaN

    fdmax   = np.ma.masked_where(np.isnan(fdmax),fdmax  )
    kdmax   = np.ma.masked_where(np.isnan(fdmax),kdmax  )
    dirmax  = np.ma.masked_where(np.isnan(fdmax),dirmax )
    kerrmax = np.ma.masked_where(np.isnan(fdmax),kerrmax)
    maskw   = fdmax.mask
    
    #Matlab like flaten for comparison
    xwcdf  = xwd.flatten(1)
    maskwf = maskw.flatten(1)
    xwcd   = xwcdf[~maskwf]

    ywcdf  = ywd.flatten(1)
    ywcd   = ywcdf[~maskwf]

    fwdf   = fdmax.flatten(1)
    fwd    = fwdf[~maskwf]
    sig    = 2 * np.pi * fwd

    k_finf   = kdmax.flatten(1)
    k_fin    = k_finf[~maskwf]
    
    skwdf   = kerrmax.flatten(1)
    skwd    = skwdf[~maskwf]
else:
    ###synthetic
    ncf_wdata=data_dir+'obs_wave/syn/wav_syn_obs.nc'
    print 'SYN data obs file name > ', ncf_wdata
    nc_wdata=netCDF4.Dataset(ncf_wdata)
    ncv_wdata=nc_wdata.variables
    xwd   = np.squeeze(ncv_wdata['x']   [::jumpp])
    ywd   = np.squeeze(ncv_wdata['y']   [::jumpp])
    uwd   = np.squeeze(ncv_wdata['velx'][0,::jumpp,::jumpp])
    vwd   = np.squeeze(ncv_wdata['vely'][0,::jumpp,::jumpp])
    dwd   = np.squeeze(ncv_wdata['depw'][0,::jumpp,::jumpp])
    Hwd   = np.squeeze(ncv_wdata['hsig'][0,::jumpp,::jumpp])
    mwd   = np.squeeze(ncv_wdata['mask'][0,::jumpp,::jumpp])
    dirwd = np.squeeze(ncv_wdata['dir'] [0,::jumpp,::jumpp])
    tm01  = np.squeeze(ncv_wdata['tm01'][0,::jumpp,::jumpp])
    xwd,ywd = np.meshgrid(xwd, ywd)
    
    # masking land point 
    #   only in case of twin experiment or may be remote sensing data
    mwd   = np.ma.masked_where(mwd==0,mwd)
    maskw = mwd.mask
    if False:
        # change obs data to vector
        xwcd   = xwd   [~maskw].flatten(1)
        ywcd   = ywd   [~maskw].flatten(1)
        uwcd   = uwd   [~maskw].flatten(1)
        vwcd   = vwd   [~maskw].flatten(1)
        dwcd    = dwd   [~maskw].flatten(1)
        Hwcd    = Hwd   [~maskw].flatten(1)
        dirwcd  = dirwd [~maskw].flatten(1)
        tm01c   = tm01  [~maskw].flatten(1)
    else:
        maskwd = maskw.flatten(1)
        xwcd   = xwd.flatten(1)
        ywcd   = ywd.flatten(1)
        uwcd   = uwd.flatten(1)
        vwcd   = vwd.flatten(1)
        dwcd    = dwd.flatten(1)
        Hwcd    = Hwd.flatten(1)
        dirwcd  = dirwd.flatten(1)
        tm01c   = tm01.flatten(1)         
        
        xwcd   = xwcd   [~maskwd]
        ywcd   = ywcd   [~maskwd]
        uwcd   = uwcd   [~maskwd]
        vwcd   = vwcd   [~maskwd]
        dwcd   = dwcd   [~maskwd]
        Hwcd   = Hwcd   [~maskwd]
        dirwcd = dirwcd [~maskwd]
        tm01c  = tm01c  [~maskwd]        
    
    fwd    = 1.0/(tm01c+0.0001)
    sig  = 2* np.pi * fwd

    print ' > Wave relative direction to current assumed to be Nautical >> Hrad Coded'
    alpha  = 270-dirwcd
    beta   = np.rad2deg(np.arctan2(vwcd,uwcd))
    gamma  = np.deg2rad(beta-alpha)  # angle btween wave and current
    kguess = approxDispersion(sigma=sig,h=dwcd)

    #Calculate wavenumber
    if False:
        #If Hwd = 0 then simple wave dispersion will be used
        #test simple dispersion
        Hwcd = np.zeros_like(Hwd)
        uwcd = uwcd * 0.0
        vwcd = vwcd * 0.0
    
    
    k_fin = Dispersion(sigma = sig,dep = dwcd, u = uwcd ,
                       v = vwcd ,ang = gamma ,kguess = kguess ,
                       Hs = Hwcd)
    #
    # assign measurement uncertainty
    skwd = np.ones(k_fin.shape) * np.random.randn(k_fin.shape[0]) * 1e-2


# Create container for each variables
ndata  = len(xwcd)
k      = component(num=ndata)
# Filling the containers
k.x    = xwcd
k.y    = ywcd
k.s    = skwd  
k.data = k_fin
k.f    = fwd

meas = {}
meas.update({'k':k})

#sys.exit()
###############################################
print '> Interpolate ensemble to obs-points'
###############################################
for field in meas.keys():
    obs=meas[field]
    print ' > Interpolation of the members for > ',field
    if field=='k':
        nobs=len(obs.x)
        hi=np.zeros((nobs,Nw),dtype='float')
        ui=np.zeros_like(hi)
        vi=np.zeros_like(hi)
        Hi=np.zeros_like(hi)
        alphai=np.zeros_like(hi)
        
        for n in range(Nw):
            hi[:,n]=interpg(xw[~maskm],yw[~maskm],dw[~maskm,n],obs.x,obs.y)
        print '  > k1/5'

        for n in range(Nw):
            ui[:,n]=interpg(xw[~maskm],yw[~maskm],uw[~maskm,n],obs.x,obs.y)
        print '  > k2/5'

        for n in range(Nw):
            vi[:,n]=interpg(xw[~maskm],yw[~maskm],vw[~maskm,n],obs.x,obs.y)
        print '  > k3/5'

        for n in range(Nw):
            Hi[:,n]=interpg(xw[~maskm],yw[~maskm],Hw[~maskm,n],obs.x,obs.y)
        print '  > k4/5'

        for n in range(Nw):
            alphai[:,n]=interpg(xw[~maskm],yw[~maskm],dirw[~maskm,n],obs.x,obs.y)
        print '  > k5/5'
        
        print ' > Wave relative direction to current assumed to be Nautical >> Hrad Coded'
        alphai  = 270-alphai
        beta  = np.rad2deg(np.arctan2(vi,ui))
        gamma = np.deg2rad(beta-alphai) # angle btwn wave and curr
        
        ui   [np.isnan(ui)]=0
        vi   [np.isnan(vi)]=0
        Hi   [np.isnan(Hi)]=0
        gamma[np.isnan(gamma)]=np.pi/2.0
        
        sigma2 = np.repeat(a=sig, repeats=Nw, axis=0).reshape(sig.shape[0],Nw)
        print '  >>> in approxDispersion() '
        kguess = approxDispersion(sigma=sigma2,h=hi)
        # discard locations for which no ensemble wavenumber estimate is
        # possible
        if True:
            hi[hi>1e3]=np.nan
            [ind]=np.where(~np.isnan(sigma2.sum(1)*hi.sum(1)*kguess.sum(1)) & ~(hi.sum(1)<=0.0) & ~(kguess.sum(1)<=0.0))
            #ind=find(~isnan(sigma.*sum(hi,2).*sum(kguess,2)) &  sum(hi<0.1,2)==0 & sum(kguess<0,2)==0);  #matlab orig ones
            tmp=obs
            obs.x    = tmp.x[ind]
            obs.y    = tmp.y[ind]
            obs.s    = tmp.s[ind]
            obs.f    = tmp.f[ind]
            obs.data = tmp.data[ind]
            sigma2   = sigma2[ind]
            hi       = hi[ind]
            ui       = ui[ind]
            vi       = vi[ind]
            kguess   = kguess[ind]
            gamma    = gamma[ind]
            Hs       = Hi[ind]
        if False:
            #test simple dispersion
            Hs = np.zeros_like(gamma)
            #Hs = Hi
            ui = ui * 0.0
            vi = vi * 0.0
        
        #[nyi,nxi]=k_all.shape    
        #obs.model=k_all.reshape(nyi*nxi)
        print '  >>> call calc_k()'
        start_calc_k     = datetime.now()
        multi            = True
        if multi:
            npoint,nmemw = ui.shape
            nbr_chunks   = 10
            chunk_size   = npoint / nbr_chunks
            chunks = [(sigma2[x*chunk_size:(x+1)*chunk_size,:], \
                       kguess[x*chunk_size:(x+1)*chunk_size,:], \
                       hi    [x*chunk_size:(x+1)*chunk_size,:], \
                       ui    [x*chunk_size:(x+1)*chunk_size,:], \
                       vi    [x*chunk_size:(x+1)*chunk_size,:], \
                       gamma [x*chunk_size:(x+1)*chunk_size,:], \
                       Hs    [x*chunk_size:(x+1)*chunk_size,:],)\
                       for x in xrange(nbr_chunks)]
            p = multiprocessing.Pool()
            po = p.map_async(calc_k, chunks)
            results = po.get()
            obs.model=np.zeros_like(ui)
            x=0
            for res in results:
                 obs.model[x*chunk_size:(x+1)*chunk_size,:]= res
                 x+=1        
        else:
            chunks=(sigma2,kguess,hi,ui,vi,gamma,Hs)
            obs.model = calc_k(chunks)
        
        end_calc_k=datetime.now()
        delta_calc_k=end_calc_k-start_calc_k
        print 'total_second= ', delta_calc_k.total_seconds()
        meas['k'] = obs

# remove any nan possibly created during interpolation
for field in meas.keys():
    obs=meas[field]
    failind   = np.array(np.where(np.isnan(obs.model.sum(1)))).squeeze()
    keepind   = np.setdiff1d(range(len(obs.x)),failind)
    obs.x     = obs.x    [keepind]
    obs.y     = obs.y    [keepind]
    obs.s     = obs.s    [keepind]
    obs.data  = obs.data [keepind]
    obs.model = obs.model[keepind,:]

##################################################
# prepare output
##################################################
for field in meas.keys():
    obs        = meas[field]
    namep      = '/wav_'+ field+'.nc'
    out_name   = wav_member_dir + namep
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
    outnc.history  = 'moghimis@gmail.com  wave data structure for assimilation '+datetime.now().isoformat()
    #outnc.att=  infile[-12:]
    outnc.close()

if False:
    import cPickle as pickle
    pickle.dump( meas, open(pick_name , "wb" ) )
    print 'Wave pickle is ready at > ', pick_name

args = sys.argv
scr_name = args[0]    
scr_dir1 = os.getcwd()
os.system('cp -fr  ' + scr_name + '    ' + wav_member_dir)
os.system('cp -fr  base_info.py        ' + wav_member_dir)

print '  END > '
