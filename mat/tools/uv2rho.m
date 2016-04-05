function out=uv2rho(in,type,dir)
%
% out=uv2rho(in,type)
%
% Maps roms u (or v) grid to rho-grid.  Assumes standard format for
% depth-averaged roms outputs: size(in)=[nt,Mm,Lm].
%
% To map u-grid to rho-grid, set type='u'.
% To map v-grid to rho-grid, set type='v'.
%
% Ghost points outside the input grid are set to NaN
%
% to map from uv to rho, set dir=-1 (default: dir=1)
%

if(exist('dir')~=1)
  dir=1;
end

% if dealing with 3D inputs, lump time and vertical coord into one
% dimension.  Extract it back out at the end
if(length(size(in))==4)
  [nt,nz,ny,nx]=size(in);
  in=reshape(in,[nt*nz ny nx]);
elseif(length(size(in))==3)
  [nt,ny,nx]=size(in);
  nz=1;
elseif(length(size(in))==2)
  nz=1;
  nt=1;
  [ny,nx]=size(in);
else
  error(['can''t work with arrays with ' num2str(length(size(in))) ' dims'])
end

% it's easier to work with data in permuted form...
if(length(size(in))==3)
  in=permute(in,[2 3 1]);
end

if(type=='u')
  if(dir==1)
    out(1:ny,2:nx,1:(nt*nz))=(1/2)*(in(:,1:nx-1,1:nt*nz)+in(:,2:nx,1:nt*nz));
    out(1:ny,1,1:(nt*nz))=nan;
    out(1:ny,nx+1,1:nt*nz)=nan;
  else
    out(1:ny,1:nx-1)=(1/2)*(in(:,1:nx-1,1:nt*nz)+in(:,2:nx,1:nt*nz));
  end
elseif(type=='v')
  if(dir==1)
    out(2:ny,1:nx,1:(nt*nz))=(1/2)*(in(1:ny-1,:,1:nt*nz)+in(2:ny,:,1:nt*nz));
    out(1,1:nx,1:(nt*nz))=nan;
    out(ny+1,1:nx,1:(nt*nz))=nan;
  else
    out(1:ny-1,1:nx,1:(nt*nz))=(1/2)*(in(1:ny-1,:,1:nt*nz)+in(2:ny,:,1:nt*nz));
  end
end

% unpermute
if(length(size(out))==3)
  out=permute(out,[3 1 2]);
end

% unwrap time and vertical coords
[ntz,ny,nx]=size(out);
if(nz>1)
  out=reshape(out,[nt nz ny nx]);
end
