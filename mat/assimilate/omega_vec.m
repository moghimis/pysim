function Cpp = omega_vec(Cpp,dist,L)
%
% Cpp = omega_vec(Cpp,dist,L)
%
% helper function: decorrelation for long length scales (matrix to be
% multiplied onto prior covariance, element by element
%

if(L==0) return; end

a = sqrt(10/3)*L;
b = dist;
c = 0*Cpp;

boa = b./a;
ind=find(0<=b&b<=a);
c(ind) = -(1/4)*boa(ind).^5 ...
         +(1/2)*boa(ind).^4 ...
         +(5/8)*boa(ind).^3 ...
         -(5/3)*boa(ind).^2 ...
         +1;

ind=find(a<b&b<=2*a);
c(ind) = (1/12)*boa(ind).^5 ...
         -(1/2)*boa(ind).^4 ...
         +(5/8)*boa(ind).^3 ...
         +(5/3)*boa(ind).^2 ...
         -5*boa(ind) ...
         +4 ...
         -(2/3)*boa(ind).^(-1);

Cpp=Cpp.*c;


