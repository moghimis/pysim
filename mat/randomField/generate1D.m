function A=generate1D(nx,dx,Lx,Lz,N)
%
% A=generate1D(nx,dx,Lx,Lz,N)
%
% Uses method of Evensen (1994) for generating random 1D fields.  The
% output is an ensemble of N realizations, having covariance
%
%     C(dx) = Lz*exp(-3*(dx/Lx)^2)
%
% This means the covariance is roughly zero at distances greater than
% |Lx,Ly|, and the covariance is equal to exp(-1) at about |Lx,Ly|/3.
%

% Development Notes:
%
% Tested against transects from generate2D.m, looks good.  Spectral
% properties are the same.
%
% Covariance of the output matches with the specified covariance function.
% e.g. use this code (again, similar to test of generate2D.m)...
%
%  >> nx=500; dx=5; lx=30; N=250;
%  >> A=generate1D(nx,dx,lx,1,N);
%  >> C=myCov(A,A);
%  >> delta=([1:100]-50)*dx;
%  >> plot(delta,C(1:100,50)), hold on
%  >> plot(delta,exp(-3*(delta/lx).^2),'r')
%

%-------------------------------------------
% input params
%-------------------------------------------

% dx=5;
% nx=100;
% N=250;  % no. of samples
% Lx=50;  % decorr. length
verbose=0;
doperiodic=0;
useSVD=0;  % resample using SVD

%-------------------------------------------
% derived params
%-------------------------------------------
if(verbose)
  disp('defining parameters')
end

% to get a nonperiodic ensemble, define extra "ghost" gridpoints
if(doperiodic)
  n1=nx;
else
  n1=round(1.2*nx);
end

% enforce even record length
if(doperiodic & mod(n1,2)~=0)
  error('can''t do periodic output with odd record lengths')
end
n1=n1+mod(n1,2);

% define constants
pi2=2.0*pi;
deltak=pi2/(n1*dx);
kappa=pi2/((n1)*dx);
kappa2=kappa^2;

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

%-------------------------------------------
% solve systems for r1,r2,c.
%-------------------------------------------

if(verbose)
  disp('solving for r1,r2,c')
end

% define wavenumber indeces p,l, excluding l==0
l=[(-n1/2+1):(n1/2)];
ind=setdiff(1:length(l(:)),find(l==0));
ln0=l(ind);

% solve for coefficients
e = @(r1) exp( -2.0*( kappa2*ln0.^2/r1^2 ) );
f = @(r1) sum( e(r1) .* ( cos(kappa *ln0*rx) - exp(-1) ) );
r = fsolve(@(r)[f(r(1))],3./rx,optimset('Display','off'));
r1=r(1);
summ=sum( exp(-2*(kappa2*l.^2/r1^2)) );
summ=summ-1.0;
c=sqrt(1.0/(deltak*summ));

if(verbose)
  disp(['pseudo2D: r1 = ' num2str(r1)]);
  disp(['pseudo2D:  c = ' num2str(c )]);
end

% define aij matrices.  Note rotation is not enabled in this code
a11=1.0/r1^2;

%-------------------------------------------
% apply inverse transform to get realizations
%-------------------------------------------

if(verbose)
  disp(['computing ' num2str(nreal) ' realizations'])
end

% define wavenumber indeces following matlab ifft2 convention
l=[0:(n1/2)];

% define amplitudes 'C', in 1st quadrant
e=exp(-( a11*kappa2*l.^2  ));
clear C
C([0:(n1/2)]+1)=e*c*sqrt(deltak);
C([0]+1)=0;

% for each wavenumber (l) of each sample (j=1..N).  Could use parfor
% here, but it actually appears to be slower for my test case
for nn=1:nreal
  qhat =zeros(n1,1);
  qhat2=zeros(n1,1);

  % 1st quadrant: phase is arbitrary
  phi=2*pi*rand(size(C));
  phi([n1/2]+1)=0;
  qhat([0:n1/2]+1)=C.*exp(sqrt(-1)*(phi));

  % 2nd quadrant: conjugate symmetry
  for i=(n1/2+2):n1
    qhat(i) = conj(qhat(mod(n1-i+1,n1)+1));
  end

  % Invert the fourier transform to get the sample
  A(:,nn)=ifft(qhat,'symmetric')*n1;

end

% cut down to desired size
A=A(1:nx,:);

%-------------------------------------------
% Resample using SVD decomposition.  Follows method of EnKF-matlab routine
% "redraw_ensemble.m", available from NERSC EnKF website (Pavel Sakov, Geir
% Evensen).
%-------------------------------------------
if(useSVD)
  error('not yet adapted for 1D code')
  
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
  A=A(1:nx,1:N);
end

%-------------------------------------------
% correct mean and variance
%-------------------------------------------

A=A-repmat(mean(A,2),[1 N]);
A=A./repmat(std(A,[],2),[1 N])*Lz;

