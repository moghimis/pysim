%-----------------------------------------------------
% helper function: compute distances matrix for cross-covariance between
% field on (x,y) and field on (xo,yo)
%-----------------------------------------------------
function dist = covDist(x,y,xo,yo)

N=length(x);
M=length(xo);
dx = repmat(x(:),[1 M]) - repmat(xo(:)',[N 1]);
clear x xo
dy = repmat(y(:),[1 M]) - repmat(yo(:)',[N 1]);
clear y yo
dist = sqrt(dx.^2+dy.^2);

% % for very large memory cases, the following script might be better, 
% % as it avoids using repmat
% N=length(x);
% M=length(xo);
% for i=1:N
%   dist(i,:)=sqrt((x(i)-xo).^2+(y(i)-yo).^2);
% end