function [z,sc,Cs]=my_scoord(fn,kgrid)
%
% [z,sc,Cs]=my_scoord(fn,kgrid,ind)
%
% wrapper for scoord.m: generates z-coordinates from ocean_his.nc file
%
% fn: ocean_his.nc filename
% kgrid: =0 --> depths of RHO-points
%        =1 --> depths of W-points
% note that u,v and density are computed on RHO-points, boundary conditions
% are applied at W-points ( that is, h+zeta = diff(z_w([1 end])) ).
%

%--------------------------------------------------------
% original scoord.m header:
%--------------------------------------------------------
%
% [z,sc,Cs]=scoord(h,theta_s,theta_b,Tcline,N,kgrid,column,index,plt)
%
% This routine computes the depths of RHO- or W-points for a grid section
% along columns (ETA-axis) or rows (XI-axis).
%
% On Input:
%
%    h         Bottom depth (m) of RHO-points (matrix).
%    theta_s   S-coordinate surface control parameter (scalar):
%                [0 < theta_s < 20].
%    theta_b   S-coordinate bottom control parameter (scalar):
%                [0 < theta_b < 1].
%    Tcline    Width (m) of surface or bottom boundary layer in which
%              higher vertical resolution is required during streching
%              (scalar).
%    N         Number of vertical levels (scalar).
%    kgrid     Depth grid type logical switch:
%                kgrid = 0   ->  depths of RHO-points.
%                kgrid = 1   ->  depths of W-points.
%    column    Grid direction logical switch:
%                column = 1  ->  column section.
%                column = 0  ->  row section.
%    index     Column or row to compute (scalar):
%                if column = 1, then   1 <= index <= Lp
%                if column = 0, then   1 <= index <= Mp
%    plt       Switch to plot scoordinate (scalar):
%                plt = 0     -> do not plot.
%                plt = 1     -> plot.
%
% On Output:
%
%    z       Depths (m) of RHO- or W-points (matrix).
%    sc      S-coordinate independent variable, [-1 < sc < 0] at
%            vertical RHO-points (vector).
%    Cs      Set of S-curves used to stretch the vertical coordinate
%            lines that follow the topography at vertical RHO-points
%              (vector).


h       = nc_varget(fn,'h');       % bottom depth (m) on RHO-points
zeta    = nc_varget(fn,'zeta');    % free surface position
theta_s = nc_varget(fn,'theta_s'); % S-coordinate surface control parameter
theta_b = nc_varget(fn,'theta_b'); % S-coordinate bottom control parameter
Tcline  = nc_varget(fn,'Tcline');  % width (m) of surface or bottom boundary
                                   % layer in which higher vertical resolution
                                   % is required during streching
[dnames,dsizes]=my_nc_dim(fn);
N = dsizes(findCellStr(dnames,'N')); % number of vertical levels
L = dsizes(findCellStr(dnames,'xi_rho'));
M = dsizes(findCellStr(dnames,'eta_rho'));

% get all columns
column=1;
plt=0;
for i=1:size(zeta,1)
  for row=1:L
    [z(i,:,row,:),sc(i,:,row),Cs(i,:,row)]=...
      scoord(h,theta_s,theta_b,Tcline,N,kgrid, ...
             column,row,plt,squeeze(zeta(i,:,:)));
  end
end

% reorder dimensions, to match those of 3d output variables
z=permute(z,[1 4 2 3]);

% % test code:
% for row=1:L
%   [z(:,row,:),sc(:,row),Cs(:,row)]=...
%     scoord(h,theta_s,theta_b,Tcline,N,kgrid, ...
%     column,row,plt);
% end
