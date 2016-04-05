function Amat=pseudo2d(nx,ny,lde,rx,ry,dx,dy,n1,n2,theta,verbose)
%
% pseudo2d(nx,ny,lde,rx,ry,dx,dy,n1,n2,theta,verbose)
%
%
% This routine calculates the pseudo random fields using the procedure
% outlined in Evensen (1994)
%

if (lde < 1)    error('pseudo2D: error lde < 1'); end
if (rx <= 0.0)  error('pseudo2D: error, rx <= 0.0'); end
if (ry <= 0.0)  error('pseudo2D: error, ry <= 0.0'); end
if (n1 < nx)    error('pseudo2D: n1 < nx'); end
if (n2 < ny)    error('pseudo2D: n2 < ny'); end

pi2=2.0*pi;
deltak=pi2^2/((n1*n2)*dx*dy);
kappa=pi2/((n1)*dx);
kappa2=kappa^2;
lambda=pi2/((n2)*dy);
lambda2=lambda^2;
scale=1.0;

%-------------------------------------------
% solve systems for r1,r2,c.
%-------------------------------------------

% define wavenumber indeces p,l
p=[(-n2/2+1):(n2/2)];
l=[(-n1/2+1):(n1/2)];
[p,l]=meshgrid(p,l);

% solve for coefficients
e = @(r1,r2) exp( -2.0*( kappa2*l.^2/r1^2 + lambda2*p.^2/r2^2 ) );
f = @(r1,r2) sum(sum( e(r1,r2) .* ( cos(kappa *l*rx) - exp(-1) ) ));
g = @(r1,r2) sum(sum( e(r1,r2) .* ( cos(lambda*p*ry) - exp(-1) ) ));
r = fsolve(@(r)[f(r(1),r(2));g(r(1),r(2))],3./[rx ry]);
r1=r(1); r2=r(2);
summ=sum(sum( exp(-2*(kappa2*l.^2/r1^2 + lambda2*p.^2/r2^2)) ));
summ=summ-1.0;
c=sqrt(1.0/(deltak*summ));

if (verbose)
  disp(['pseudo2D: r1=  ' num2str(r1)]);
  disp(['pseudo2D: r2=  ' num2str(r2)]);
  disp(['pseudo2D:  c=  ' num2str(c )]);
end

% Rotation to angle theta
a11tmp=1.0/r1^2;
a22tmp=1.0/r2^2;
torad=-pi/180.0;
a11=a11tmp*cos(theta*torad).^2 + a22tmp*sin(theta*torad).^2;
a22=a11tmp*sin(theta*torad).^2 + a22tmp*cos(theta*torad).^2;
a12=(a22tmp-a11tmp).*cos(theta*torad).*sin(theta*torad);

% for each wavenumber (p,l) of each sample (j=1..lde)...
p=0:(n2-1);
l=0:(n1-1);
[p,l]=meshgrid(p,l);
for j=1:lde

  % Calculating the random wave phases
  phi=pi2*rand(size(p));

  % Calculating the complex wave amplitudes 'x'
  e=exp(-( a11*kappa2*l.^2 + 2.0*a12*kappa*lambda*l.*p + a22*lambda2*p.^2 ));
  fampl(:,:,1)=e.*cos(phi)*sqrt(deltak)*c;
  fampl(:,:,2)=e.*sin(phi)*sqrt(deltak)*c;
  fampl(1,1,:)=0.0;
  x=fampl(:,:,1)+sqrt(-1)*fampl(:,:,2);

  % Invert the fourier transform to get the sample
  Amat(:,:,j)=fft2(fftshift2(x),'symmetric');

end

