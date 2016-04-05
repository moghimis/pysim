function k=approxDispersion(sigma,h)
%
% k=approxDispersion(sigma,h)
%
% Approximation due to Guo (2002), see discussion here:
%  http://www.johndfenton.com/Papers/Dispersion-Relation.pdf
%
% Inputs should be of like dimension
%

g=9.8126;
k = sigma.^2/g.*(1-exp(-(sigma.*sqrt(h/g)).^(5/2))).^(-2/5);

