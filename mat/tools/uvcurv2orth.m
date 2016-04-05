function [uout,vout,alph,beta] = uvcurv2orth(x,y,uin,vin,direction)
%
% [uout,vout] = uvcurv2orth(x,y,uin,vin,[direction])
%
% Converts 2d vector directions from curvilinear grid directions to
% orthogonal directions on the curvilinear grid points
%
% If there are multiple snapshots to be converted, make input in the form
% size(uin)=[nt,ny,nx] (standard form for roms output).  Similar for vin.
%
% Optionally, set direction=-1 to go from orthogonal to curvilinear
% (i.e. the "inverse" transform); default is direction=1 ("forward").
%

if(~exist('direction')) direction=1; end

% easier to work with data in permuted form
if(length(size(uin))==3)
  uin=permute(uin,[2 3 1]);
  vin=permute(vin,[2 3 1]);
end
[ny,nx,nt]=size(uin);

% calculate grid spacing in curvilinear coords 
dxvert = diff(x,1,1);
dyvert = diff(y,1,1);
dxhorz = diff(x,1,2);
dyhorz = diff(y,1,2);

% assume constant spacing at mesh boundary
dxvert = cat(1,dxvert(1,:),dxvert);
dyvert = cat(1,dyvert(1,:),dyvert);
dxhorz = cat(2,dxhorz(:,1),dxhorz);
dyhorz = cat(2,dyhorz(:,1),dyhorz);

% compute local angles between curvilinear and orthogonal x and y angles
% (alph and beta, respectively)
alph=atan2(dyhorz,dxhorz);
beta=atan2(dxvert,dyvert);

% apply forward transform
if(direction==1)

  a=cos(alph);
  b=sin(beta);
  c=sin(alph);
  d=cos(beta);

  % handle multiple snapshots
  a=repmat(a,[1 1 nt]);
  b=repmat(b,[1 1 nt]);
  c=repmat(c,[1 1 nt]);
  d=repmat(d,[1 1 nt]);

  uout = a.*uin + b.*vin;
  vout = c.*uin + d.*vin;

% apply inverse transform
else

  %   % test code: direct numerical inverse of forward transform
  %   for i=1:size(uin,1)
  %     for j=1:size(uin,2)
  %       b = [uin(i,j) vin(i,j)]';
  %       A = [[cos(alph(i,j)) sin(beta(i,j))]
  %            [sin(alph(i,j)) cos(beta(i,j))]];
  %       % note b=A*c is the forward transform
  %       c=A\b;
  %       uout(i,j)=c(1);
  %       vout(i,j)=c(2);
  %     end
  %   end


  % handle multiple snapshots
  alph=repmat(alph,[1 1 nt]);
  beta=repmat(beta,[1 1 nt]);
  
  % regular code: analytical inverse
  denom=-cos(alph).*cos(beta)+sin(alph).*sin(beta);

  uout =-(uin.*cos(beta)-vin.*sin(beta))./denom;
  vout = (uin.*sin(alph)-vin.*cos(alph))./denom;

end

% unpermute
if(length(size(uin))==3)
  uout=permute(uout,[3 1 2]);
  vout=permute(vout,[3 1 2]);
end
