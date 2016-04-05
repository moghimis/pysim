from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pre-process ROMS and velocity data and put them in vector form.
 - A little bit of cleaning as well
 - Finally interpolated data and model will be save in a netcdf file to read by assimilation routine 
#
#
    NEED TO CHECK:
    > 1. Division by numbers 4/5 check not to devide 
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
jumpp     = base_info.jump_cur
real_curr = base_info.real_data

#pysim_inp = base_dir + '/pysim_inp.txt'
####################################################
#Set input parameters
#fdata1 = open(pysim_inp)
#for  line in fdata1.readlines():
#    print line
#    if 'itr'    in line: itr          = int(line.split()[-1]) 
#fdata1.close()
#####################################################
inp_dir  = base_dir + '/inp/'
run_id   = '/run_'+str(1000+itr) 
#### Funcs
##   interpolation
methodi='csa'
if methodi=='csa':
    import octant.csa as csa
elif methodi=='tri':
    from delaunay import triangulate

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
    
class component:
    def __init__(self,num):
        self.x=np.zeros((num),dtype='float')
        self.y=np.zeros((num),dtype='float')
        self.data=np.zeros((num),dtype='float')
        self.s=np.zeros((num),dtype='float')
#### End of Funcs
cur_member_dir = base_dir+run_id+'/04_mem_adj'
data_dir       = inp_dir+'/'
#### Pickle name and check if exist !
namep      = '/cur_jump-'+str(jumpp)+'__itr-'+str(itr)+'.p'
pick_name  = cur_member_dir + namep
#if  os.path.isfile(pick_name):
#    sys.exit('CUR file exist !')
###########################################################
# Read circulation  data
# load the data grid for currents 
#ncf_cdata=data_dir+'obs/obs_every2ocean_avg_region.nc'
if real_curr:
    ncf_cdata = data_dir+'obs/sar/uASAR.nc'
    print "Real cur data > ", ncf_cdata
    nc_cdata  = netCDF4.Dataset(ncf_cdata)
    ncv_cdata = nc_cdata.variables
    xcd  = np.squeeze(ncv_cdata['x'][:])
    ycd  = np.squeeze(ncv_cdata['y'][:])
    ucd  = np.squeeze(ncv_cdata['u'][0,:])
    vcd  = np.squeeze(ncv_cdata['v'][0,:])
    sucd = 1.0 * np.squeeze(ncv_cdata['u_err'][:])
    svcd = 1.0 * np.squeeze(ncv_cdata['v_err'][:])
    #zcd=np.squeeze(ncv_cdata['z'][0,:])
    #mcd=np.squeeze(ncv_cdata['mask'][:])
    maskd = ucd > 1e5 
    ucd   = np.ma.masked_array(ucd,maskd)
    xcd   = np.ma.masked_array(xcd,maskd)
    ycd   = np.ma.masked_array(ycd,maskd)
    ucd   = np.ma.masked_array(ucd,maskd)
    vcd   = np.ma.masked_array(vcd,maskd)
    sucd  = np.ma.masked_array(sucd,maskd)    
    svcd  = np.ma.masked_array(svcd,maskd)        
    
    xcdm  = xcd.flatten(1)
    ycdm  = ycd.flatten(1)
    ucdm  = ucd.flatten(1)
    vcdm  = vcd.flatten(1)
    sucdm = sucd.flatten(1)
    svcdm = svcd.flatten(1)
    
    xcdm = xcdm.compressed()  
    ycdm = ycdm.compressed()
    ucdm = ucdm.compressed()
    vcdm = vcdm.compressed() 
    sudm = sucdm.compressed()
    svdm = svcdm.compressed()
else:
    ncf_cdata=data_dir+'obs/syn/syn1nri_his.nc'
    print "SYN cur data > ", ncf_cdata
    nc_cdata=netCDF4.Dataset(ncf_cdata)
    ncv_cdata=nc_cdata.variables
    xcd = np.squeeze(ncv_cdata['x_rho']  [:])
    ycd = np.squeeze(ncv_cdata['y_rho']  [:])
    ucd = np.squeeze(ncv_cdata['u']    [0,:])
    vcd = np.squeeze(ncv_cdata['v']    [0,:])
    maskcd = np.squeeze(ncv_cdata['mask'][:])
    maskd = (maskcd==0)
    xcd,ycd = np.meshgrid(xcd, ycd)
    sucd = 0.04 * np.ones(ucd.shape) * np.random.randn(ucd.shape[0],ucd.shape[1])
    svcd = 0.04 * np.ones(vcd.shape) * np.random.randn(vcd.shape[0],vcd.shape[1])
    #Matlab like fllaten for comparison
    xcdf= xcd.flatten(1)
     
    maskdf = maskd.flatten(1)
    xcdm = xcdf[~maskdf]
     
    ycdf= ycd.flatten(1)
    ycdm = ycdf[~maskdf]
     
    ucdf= ucd.flatten(1)
    ucdm = ucdf[~maskdf]
     
    vcdf= vcd.flatten(1)
    vcdm = vcdf[~maskdf]
     
    sucdf = sucd.flatten(1)
    sudm   = sucdf[~maskdf]
     
    svcdf = svcd.flatten(1)
    svdm   = svcdf[~maskdf]
    
#     # masking land point 
if jumpp >1:
    ucdm    = ucdm[::jumpp]
    vcdm    = vcdm[::jumpp]
    xcdm    = xcdm[::jumpp]
    ycdm    = ycdm[::jumpp]
    sudm    = sudm[::jumpp]
    svdm    = svdm[::jumpp]

# Create container for each variables
ndata=len(xcd)
u=component(num=ndata)
v=component(num=ndata)
#z=component(num=ndata)

# Filling the containers
u.x    = xcdm
u.y    = ycdm
u.s    = sudm
u.data = ucdm

v.x    = xcdm
v.y    = ycdm
v.s    = svdm
v.data = vcdm

#z.x=xcd
#z.y=ycd
#z.s=szcd
#z.data=zcd

meas={}
meas['u']=u
meas['v']=v

######################
print 'Loading ROMS members '
######################
memdir  = cur_member_dir+'/member1*'
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
xc    = ncv['x_rho'][:];xc1=xc[1,:]
yc    = ncv['y_rho'][:];yc1=yc[:,1]
[ny,nx]=h.shape
nc.close()

uc    = np.zeros((ny,nx,nmemc),dtype='float')
vc    = np.zeros_like(uc)
zc    = np.zeros_like(uc)
hc    = np.zeros_like(uc)
mc    = np.zeros_like(uc)

failind=[]
for i in range(nmemc):
    if np.mod(i,50)==0: print ' > CUR member read in ', (i*1.0/nmemc)*100.0, '%'
    try:
        #ncf=dirlist[i]+'/nri_his1000.nc'
        ncf   = glob.glob(dirlist[i]+'/nri_h*')[0]

        nc=netCDF4.Dataset(ncf)
        ncv=nc.variables
        hc[:,:,i]=ncv['h'][:,:]
        uc[:,:,i]=ncv['ubar'][:,:]
        vc[:,:,i]=ncv['vbar'][:,:]
        zc[:,:,i]=ncv['zeta'][:,:]
        mc[:,:,i]=ncv['mask'][:,:]
        nc.close()
    except IOError:
        failind.append(i)

# In case some members failed for whatever reason: discard thempcol    
if(len(failind) > 0):
  print 'WARNING WARNING WARNING: the following members failed:'
  for i in failind:
    print 'member -> ', dirlist[i]
  
  keepind=np.setdiff1d( range(nmemc),failind)
  uc=uc[:,:,keepind]
  vc=vc[:,:,keepind]
  zc=zc[:,:,keepind]
  mc=mc[:,:,keepind]
  hc=hc[:,:,keepind]
  dirlistc=dirlist[keepind]

#mask out land points
temp = uc.sum(2) + vc.sum(2) + hc.sum(2) +zc.sum(2)
uc = np.ma.masked_where(mc==0 | np.isnan(temp) , uc)
vc = np.ma.masked_where(mc==0 | np.isnan(temp) , vc)
hc = np.ma.masked_where(mc==0 | np.isnan(temp) , hc)
zc = np.ma.masked_where(mc==0 | np.isnan(temp) , zc)
maskm = mc[:,:,0]==0
[nyc,nxc,Nc]=hc.shape

print '  > Number of CUR members >', Nc
###############################################
print 'Interpolate ensemble to obs-points'
###############################################
for field in meas.keys():
    obs=meas[field]
    print ' > Interpolation of the members for > ',field
    if field in ['u','v','z']  :
        if field=='u' : data=uc
        if field=='v' : data=vc
        if field=='z' : data=zc
        
        nobs=len(obs.x)
        datai=np.zeros((nobs,Nc),dtype='float')
        
        #xy      = np.array(zip(xc[~maskm],yc[~maskm]))
        #tmp,ind = np.unique(xy, return_index=True)
        #xcu     = xc[~maskm][ind]
        #ycu     = yc[~maskm][ind]
        
        for n in range(Nc):
            datai[:,n]=interpg(xc[~maskm],yc[~maskm],data[~maskm,n],obs.x,obs.y)
        obs.model=datai
        meas[field]=obs

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
###################################################
#adding a dummy var to make it similar to wave dictionary
for field in ['u','v']: #,'z']:
    obs=meas[field]
    obs.f=np.zeros_like(obs.x)
    meas[field]=obs
##################################################
# prepare output
##################################################
for field in ['u','v']: #,'z']:
    obs        = meas[field]
    namep      = '/cur_'+ field+'.nc'
    out_name   = cur_member_dir + namep
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

# Member prior info out
namep      = '/cur_members_prior.nc'
out_name   = cur_member_dir + namep
outnc      = netCDF4.Dataset(out_name,'w',format='NETCDF3_CLASSIC')
ny,nx,nmem = hc.shape

outnc.createDimension('ny'   , ny )
outnc.createDimension('nx'   , nx )
outnc.createDimension('nmem' , nmem )

p0 = outnc.createVariable('x_rho','f8',('ny','nx',))
p0.missing_value = -9999.0
p0[:]            = xc

p1 = outnc.createVariable('y_rho','f8',('ny','nx',))
p1.missing_value = -9999.0
p1[:]            = yc

p2 = outnc.createVariable('h_mems','f8',('ny','nx','nmem',))
p2.missing_value = -9999.0
p2[:]            = hc

p3 = outnc.createVariable('u_mems','f8',('ny','nx','nmem',))
p3.missing_value = -9999.0
p3[:]            = uc

p4 = outnc.createVariable('v_mems','f8',('ny','nx','nmem',))
p4.missing_value = -9999.0
p4[:]            = vc

outnc.history  = 'moghimis@gmail.com, mems prior for assimilation '+datetime.now().isoformat()
#outnc.att=  infile[-12:]
outnc.close()

if False:
    import cPickle as pickle
    ###
    namep      = '/cur_jump-'+str(jumpp)+'__itr-'+str(itr)+'.p'
    pick_name = cur_member_dir + namep
    pickle.dump( meas, open(pick_name , "wb" ) )
    print 'Circulation pickle is ready at > ', pick_name
    ###
    
    out_pick={'xc':xc,'yc':yc,'hc':hc}
    
    namep      = '/hp_jump-'+str(jumpp)+'__itr-'+str(itr)+'.p'
    pick_name = cur_member_dir + namep
    pickle.dump( out_pick, open(pick_name , "wb" ) )
    
args = sys.argv
scr_name = args[0]    
scr_dir1 = os.getcwd()
os.system('cp -fr  ' + scr_name + '    ' + cur_member_dir)
os.system('cp -fr  base_info.py        ' + cur_member_dir)

print 'Members CUR data is ready at > ', out_name
## end
##################################################

