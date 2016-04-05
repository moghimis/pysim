function Q=randrot(nrens)
%
% Q=randrot(nrens)
%
% This routine generates a real orthogonal random matrix.
% The algorithm is the one by Francesco Mezzadri (2007), How to generate
% random matrices from the classical compact groups, Notices of the AMS,
% Vol. 54, pp 592-604.
%
% 1. First a matrix with independent random normal numbers are simulated.
%
% 2. Then the QR decomposition is computed, and Q will then be a random
% orthogonal matrix.
%
% 3. The diagonal elements of R are extracted and we construct the
% diagonal matrix X(j,j)=R(j,j)/|R(j,j)|
%
% 4. An updated Q'=Q X is computed, and this is now a random orthogonal
% matrix with a Haar measure.
%

% construct random matrix
A=rand(nrens,nrens);
B=rand(nrens,nrens);
tiny=1e-10;
Q = sqrt(-2.*log(A+tiny)) * cos(2.*pi*B);

% perform QR decomposition and normalization
[Q,R] = qr(Q);
X=diag(diag(R));
Q=Q*X;
