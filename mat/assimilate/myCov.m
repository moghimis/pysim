function Cab = myCov(a,b)
%
% Cab = myCov(a,b)
%
% Computes sample covariance 'a' and 'b'.  Inputs should have dimensions
% as follows...
%
%    size(a) = [ M, N ]
%    size(b) = [ P, N ]
%
% That is, each matrix consists of N samples of an Mx1 (or Px1) vector.
% Note M and P need not be equal (unlike the builtin matlab function cov.m).
% Output 'Cab' is the sample covariance, which has dimensions MxP.
%

if(size(a,2) ~= size(b,2))
  error('inputs must have same number of samples (2nd dimension)')
end

da=(a-repmat(nanmean(a,2),1,size(a,2)));
db=(b-repmat(nanmean(b,2),1,size(b,2)));
N=size(a,2);
Cab = (1/(N-1))*da*db';
