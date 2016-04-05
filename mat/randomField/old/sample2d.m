function sample2d(nx,ny,nrens,nre,dx,dy,rx,ry,theta,samp_fix,verbose)
%
% This routine samples pseudo random fields with improved independence or
% orthogonality.  This is done by first drawing a large sample and construct
% the final sample using the dominant singular vectors of the large sample.
%
% Ported from fortran code by Evensen
%

error('code not finished being ported')

n1=round(nx*1.2);
n2=round(ny*1.2);
if (verbose) then
  disp(['nx=' num2str(nx)]);
  disp(['ny=' num2str(ny)]);
  disp(['n1=' num2str(n1)]);
  disp(['n2=' num2str(n2)]);
end

n=nx*ny;
ns=nre*nrens;
msx=min(ns,n);
nsx=min(nrens,n);

% Start with oversized ensemble of ns members
disp(['sample2d: calling pseudo2d'])
A=pseudo2d(nx,ny,ns,rx,ry,dx,dy,n1,n2,theta,verbose);
disp(['sample2d: pseudo2d done'])

% generate a random orthogonal matrix
VT1=randrot(nsx);

% Compute SVD of oversized ensemble
for i=1:nsx
  [UU(:,:,i),sig(:,:,i),VT(:,:,i)] = svd(A(:,:,i));
end

% Generate first min(nrens,n) members
ishape3=(/nx,ny,nsx/)
UU=reshape(U(:,1:nsx),ishape3)
A2=zeros([size(A,1) size(A,2) nsx]);
for j=1:nsx
  for i=1:nsx
    A2(:,:,j)=A2(:,:,j)+UU(:,:,i)*sig(i)/sqrt(nre)*VT1(i,j);
  end
end
disp('sample2d: improved ensemble done')

% subtract mean and correct variance
warning('todo: recenter and scale the ensemble')
