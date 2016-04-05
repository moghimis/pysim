% will be python very soon >>>
%#!/usr/bin/env python
%# -*- coding: utf-8 -*-
%"""
% Generate bathymetry ensemble
%"""
%__author__    = "Saeed Moghimi"
%__copyright__ = "Copyright 2015, Oregon State University"
%__license__   = "GPL"
%__version__   = "0.1"
%__email__     = "moghimis@gmail.com"
%
%#############################################
%# Saeed Moghimi; moghimis@gmail.com   
%# Logs:
%# 1.0 Greg Wilson adapted for matlab
%# 2.0 
%# 3.0 
%# 4.0 
%#########################################################


%
% generate bathymetry ensemble
%
%inp_dir=getenv('INP_DIR');
%addpath ([inp_dir '/scr/mat/tools/'])
%addpath ([inp_dir '/scr/mat/randomField/'])
%addpath ([inp_dir '/scr/mat/assimilate/'])

%scr_dir=getenv('SCR_DIR');

% Read params from text file
%param=load('input_param.txt')
%N    =param(1)
% perturbation scales
%Li   =param(2)
%Lj   =param(3)
%Lz   =param(4)
%type1=param(5)



fid=fopen('input_param.txt');
param=textscan(fid,'%s');
scr_dir      =        cell2mat(param{1}(1))
inp_dir      =        cell2mat(param{1}(2))
N            =str2num(cell2mat(param{1}(3)))
Li           =str2num(cell2mat(param{1}(4)))
Lj           =str2num(cell2mat(param{1}(5)))
Lz           =str2num(cell2mat(param{1}(6)))
type1        =str2num(cell2mat(param{1}(7)))
fclose(fid);

dh_min = -3
disp('dh min set to -3')

addpath ([scr_dir '/mat/tools/'])
addpath ([scr_dir '/mat/randomField/'])
addpath ([scr_dir '/mat/assimilate/'])

out_bathy='../01_bat_inp/bath'
infile='prior.nc'


if(matlabpool('size')==0)
 matlabpool
end

if(type1==1)
    disp('load prior bathymetry')
    h=nc_varget(infile,'h');
    x=nc_varget(infile,'x_rho'); x=x(1,:)';
    y=nc_varget(infile,'y_rho'); y=y(:,1);

    dx=x(2)-x(1);
    dy=y(2)-y(1);
    nx=length(x);
    ny=length(y);

    dh = generate2Dp(nx,ny,dx,dy,Li,Lj,Lz,N);

    % add the perturbations to the prior
    h0=h;
    h=repmat(h0,[1 1 N])+dh;

    [x,y]=meshgrid(x,y);
    for n=1:N
        disp(n) 
        mat2roms(x,y,dh(:,:,n),[out_bathy  num2str(1000+n) '.nc'],0);
    end
else
    grdfile='nri_curv_grid.nc'
    %N=200;

    disp('load prior bathymetry')
    h=nc_varget(infile,'h');
    x=nc_varget(infile,'x_rho'); x=x(1,:)';
    y=nc_varget(infile,'y_rho'); y=y(:,1);
    dx=x(2)-x(1);
    dy=y(2)-y(1);
    nx=length(x);
    ny=length(y);

    disp('load curvilinear grid')
    xc=nc_varget(grdfile,'x_rho');
    yc=nc_varget(grdfile,'y_rho');
    [nyc,nxc]=size(xc);
    [ic,jc]=meshgrid(1:nxc,1:nyc);

    disp('perturbation scales')
    %Li=25;      %Redius equal to LixLj times number of grid point perturbations
    %Lj=25;
    %Lz=0.5;

    % all horizontal units are in grid points
    dhc=generate2Dp(nxc,nyc,1,1,Li,Lj,Lz,N);

    disp('Interpolation on the actual grid')
    [xg,yg]=meshgrid(x,y);
    parfor n=1:N
       F=TriScatteredInterp(xc(:),yc(:),reshape(dhc(:,:,n),[nyc*nxc 1]));
       dh(:,:,n)=F(xg,yg);
    end
   
    [x,y]=meshgrid(x,y);
    for n=1:N
        disp(n) 
        mat2roms(x,y,dh(:,:,n),[out_bathy  num2str(1000+n) '.nc'],dh_min);
    end
end

 matlabpool close
