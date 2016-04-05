function [dnames,dsizes]=my_nc_dim(fname);
%
% [dnames,dsizes]=my_nc_dim(fname);
%
% GW: altered version of nc_dim which returns a cell-array of names
% rather than a single string of names
%


%-----------------------------------------------------------------------
%  Open NetCDF file.
%-----------------------------------------------------------------------

[ncid]=mexcdf('ncopen',fname,'nc_nowrite');
if (ncid == -1),
  error(['NC_DIM: ncopen - unable to open file: ', fname]);
  return
end
 
%-----------------------------------------------------------------------
%  Supress all error messages from NetCDF.
%-----------------------------------------------------------------------
 
[ncopts]=mexcdf('setopts',0);

%-----------------------------------------------------------------------
% Inquire about contents.
%-----------------------------------------------------------------------

[ndims,nvars,natts,recdim,status]=mexcdf('ncinquire',ncid);
if (status == -1),
  error(['NC_DIM: ncinquire - cannot inquire file: ',fname]);
end,

%-----------------------------------------------------------------------
% Inquire about dimensions
%-----------------------------------------------------------------------

for n=1:ndims;
  [name,size,status]=mexcdf('ncdiminq',ncid,n-1);
  if (status == -1),
    error(['NC_DIM: ncdiminq - unable to inquire about dimension ID: ',...
          num2str(n)]);
  else,
    lstr=length(name);
    dnames{n}=name(1:lstr);
    dsizes(n)=size;
  end,
end,

%-----------------------------------------------------------------------
% Close NetCDF file.
%-----------------------------------------------------------------------

[status]=mexcdf('ncclose',ncid);
if (status == -1),
  error(['NC_DIM: ncclose - unable to close file: ', fname]);
  return
end,

return