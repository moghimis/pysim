%
% test of algorithm for pseudo2d.m.  This is a precursor to generate.m,
% used for development and testing.
%
clear

%-------------------------------------------
% input params
%-------------------------------------------

dx=5;
dy=5;
nx=100;
ny=200;
N=250;  % no. of samples
Lx=50;  % decorr. length
Ly=100;
verbose=1;
doperiodic=0;
useSVD=0;  % resample using SVD
doTest=0;  % test against evensen's code

%-------------------------------------------
% derived params
%-------------------------------------------
if(verbose)
  disp('defining parameters')
end

% check for even record length
if(mod(nx,2)~=0 | mod(ny,2)~=0)
  error('this code hasn''t been extended for odd record lengths')
end

% to get a nonperiodic ensemble, define extra "ghost" gridpoints
if(doperiodic)
  n1=nx;
  n2=ny;
else
  n1=round(1.2*nx);
  n2=round(1.2*ny);
  n1=n1+mod(n1,2);
  n2=n2+mod(n2,2);
end

% define constants
pi2=2.0*pi;
deltak=pi2^2/((n1*n2)*dx*dy);
kappa=pi2/((n1)*dx);
kappa2=kappa^2;
lambda=pi2/((n2)*dy);
lambda2=lambda^2;

if(useSVD)
  beta=3;
  nreal=beta*N;  % oversized ensemble
else
  nreal=N;
end

% rescale decorrelation lengths such that we will get the following form
% for the covariance as a function of distance delta:
%
%     C(delta)=exp(-3*(delta/Lx)^2)
%
rx=Lx/sqrt(3);
ry=Ly/sqrt(3);

%-------------------------------------------
% solve systems for r1,r2,c.
%-------------------------------------------

if(verbose)
  disp('solving for r1,r2,c')
end

% define wavenumber indeces p,l, excluding p==l==0
p=[(-n2/2+1):(n2/2)];
l=[(-n1/2+1):(n1/2)];
[p,l]=meshgrid(p,l);
ind=setdiff(1:length(p(:)),find(p==0 & l==0));
pn0=p(ind);
ln0=l(ind);

% solve for coefficients
e = @(r1,r2) exp( -2.0*( kappa2*ln0.^2/r1^2 + lambda2*pn0.^2/r2^2 ) );
f = @(r1,r2) sum( e(r1,r2) .* ( cos(kappa *ln0*rx) - exp(-1) ) );
g = @(r1,r2) sum( e(r1,r2) .* ( cos(lambda*pn0*ry) - exp(-1) ) );
r = fsolve(@(r)[f(r(1),r(2));g(r(1),r(2))],3./[rx ry],...
           optimset('Display','off'));
r1=r(1); r2=r(2);
summ=sum(sum( exp(-2*(kappa2*l.^2/r1^2 + lambda2*p.^2/r2^2)) ));
summ=summ-1.0;
c=sqrt(1.0/(deltak*summ));

if(verbose)
  disp(['pseudo2D: r1 = ' num2str(r1)]);
  disp(['pseudo2D: r2 = ' num2str(r2)]);
  disp(['pseudo2D:  c = ' num2str(c )]);
end

% define aij matrices.  Note rotation is not enabled in this code
a11=1.0/r1^2;
a22=1.0/r2^2;
a12=0*a11;

%-------------------------------------------
% apply inverse transform to get realizations
%-------------------------------------------

if(verbose)
  disp(['computing ' num2str(nreal) ' realizations'])
end

% define wavenumber indeces following matlab ifft2 convention
l=[0:(n1/2)];
p=[0:(n2/2)];
[p,l]=meshgrid(p,l);

% define amplitudes 'C', in 1st quadrant
e=exp(-( a11*kappa2*l.^2 + 2.0*a12*kappa*lambda*l.*p + a22*lambda2*p.^2 ));
clear C
C([0:(n1/2)]+1,[0:(n2/2)]+1)=e*c*sqrt(deltak);
C([0]+1,:)=0;
C(:,[0]+1)=0;

% for each wavenumber (p,l) of each sample (j=1..N)...
for nn=1:nreal
  qhat =zeros(n1,n2);
  qhat2=zeros(n1,n2);

  % 1st quadrant: phase is arbitrary
  phi=2*pi*rand(size(C));
  phi(:,[n2/2]+1)=0;
  phi([n1/2]+1,:)=0;
  qhat([0:n1/2]+1,[0:n2/2]+1)=C.*exp(sqrt(-1)*(phi));

  % 3rd quadrant: phase is also arbitrary
  phi2=2*pi*rand(size(C));
  phi2(:,[n2/2]+1)=0;
  phi2([n1/2]+1,:)=0;
  qhat2([0:n1/2]+1,[0:n2/2]+1)=C.*exp(sqrt(-1)*(phi2));
  for j=(n1/2+1):(n1-1)
    for i=1:(n2/2)
      qhat(j+1,i+1)=conj(qhat2([n1-j]+1,i+1));
    end
  end
  qhat([(n1/2+1):(n1-1)]+1,1)=0;

  % 2nd and 4th quadrants are set by conjugate symmetry
  for i=(n2/2+2):n2
    for j=1:n1
      qhat(j,i) = conj(qhat(mod(n1-j+1,n1)+1,mod(n2-i+1,n2)+1));
    end
  end

  % Invert the fourier transform to get the sample
  A(:,:,nn)=ifft2(qhat,'symmetric')*n1*n2;

end

% cut down to desired size
A=A(1:nx,1:ny,:);

%-------------------------------------------
% Resample using SVD decomposition.  Follows method of EnKF-matlab routine
% "redraw_ensemble.m", available from NERSC EnKF website (Pavel Sakov, Geir
% Evensen).
%-------------------------------------------
if(useSVD)

  if(verbose)
    disp(['resampling for ' num2str(N) ' realizations'])
  end

  % % generate a random orthogonal matrix
  % VT1=genR(N);
  % 
  % % Compute SVD of oversized ensemble
  % Avec=reshape(A,[nx*ny nreal]);
  % [UU,sig,VT] = svd(Avec,0);
  % 
  % % rescale singular values
  % sig=diag(diag(sig(1:N,1:N)))*sqrt((N-1)/(nreal-1));
  % 
  % % Generate new members
  % A2vec=UU(:,1:N)*sig*VT1;
  % A2=reshape(A2vec,[nx ny N]);

  Avec=reshape(A,[nx*ny nreal]);
  B=[ones(N-1)/(sqrt(N)*(sqrt(N)+1))-eye(N-1),ones(N-1,1)/sqrt(N)];
  [Q,~]=qr(randn(N-1));
  VT=[Q*ones(N-1)/(sqrt(N)*(sqrt(N)+1))-Q,Q*ones(N-1,1)/sqrt(N);
      ones(1,N)/sqrt(N)];
  [U,L]=svd(Avec,0);
  mnmin=min([N,nx*ny,nreal]);
  if mnmin==N
    L(mnmin,mnmin)=0; % to ensure that Anew iz zero-centered
  end
  L=diag(diag(L)*sqrt((N-1)/(nreal-1)));
  A2=U(:,1:mnmin)*[L(1:mnmin,1:mnmin),zeros(mnmin,N-mnmin)]*VT;
  A2=reshape(A2,[nx ny N]);

else
  A=A(1:nx,1:ny,1:N);
end

%-------------------------------------------
% repeat using evensen's code and compare
%-------------------------------------------
if(doTest)
warning off
addpath ~/work/roms/User/tools/perturber/evensenCode/ % readOutput.m
warning on

if(verbose)
  disp('repeating using Evensen''s code')
end

% note the version 'generate2' doesn't apply recentering and variance
% corrections.  It also writes output in double precision
cmd=['/home/server/student/homes/gwilson/work/roms/User/tools/' ...
     'perturber/evensenCode/generate2'];

% execute the code
fout=['tmpout'];
args = [num2str(ny  ) ' ' ...
        num2str(nx  ) ' ' ...
        num2str(dy  ) ' ' ...
        num2str(dx  ) ' ' ...
        num2str(N ) ' ' ...
        num2str(ry  ) ' ' ...
        num2str(rx  ) ' ' ...
        num2str(0   ) ' ' ...
        num2str(fout)];
[a,b] = unix([cmd ' ' args])
Ae=readOutput(fout,ny,nx,dy,dx,N,'l','double');
delete(fout);

% compare stats from Evensen's code vs my code
subplot(2,2,1)
pcolor(std(Ae,[],3)),sf
caxis(1+[-1 1]*.2), colorbar
title('Evensen std')
subplot(2,2,3)
pcolor(mean(Ae,3)),sf
caxis([-1 1]*.25), colorbar
title('Evensen mean')
subplot(2,2,2)
pcolor(std(A,[],3)),sf
caxis(1+[-1 1]*.2), colorbar
title('My std')
subplot(2,2,4)
pcolor(mean(A,3)),sf
caxis([-1 1]*.25), colorbar
title('My mean')

% result: both codes have a problem getting convergence in the mean and
% stdev of the ensemble.  The fix in both cases is to apply recentering
% and scaling, done next

end

%-------------------------------------------
% correct mean and variance
%-------------------------------------------

A=A-repmat(mean(A,3),[1 1 N]);
A=A./repmat(std(A,[],3),[1 1 N]);

%-------------------------------------------
% reshape the output to matlab conventions
%-------------------------------------------

A=permute(A,[2 1 3]);


%-------------------------------------------
% done
%-------------------------------------------

% plot an example member
pcolor([1:nx]*dx,[1:ny]*dy,A(:,:,1)),sf

