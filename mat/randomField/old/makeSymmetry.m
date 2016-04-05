function x=makeSymmetry(x)
%
% Utility to generate a conjugate-symmetric matrix from a 2d input matrix
% x.  Used for fft2
%

[n1,n2]=size(x);
if(mod(n1,2)==0)  % even
  x=[x;flipud(conj(x))];
else
  x=[x;flipud(conj(x(1:n1-1,:)))];
end
if(mod(n2,2)==0)  % even
  x=[x fliplr(conj(x))];
else
  x=[x fliplr(conj(x(:,1:n2-1)))];
end

  