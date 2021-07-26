function [a] = mon2chebsym(B,k,j)
%CHEB2MON  Monomial to Chebyshev basis conversion.
%   A = MON2CHEB(B, k, j) converts polynomial B (of degree k) given in monomial basis to 
%   Chebyshev basis A. The polynomial must be given with its coefficients
%   in descending order, i.e. B = B_N*x^N + ... + B_1*x + B_0
%
%   Example: 
%    Suppose we have a polynomial in the monomial basis: 
%    b2*x^2 + b1*x + b0, 
%    with b2=2, b1=0, b0=-2 for example.
%    We want to express the polynomial in the Chebyshev basis
%    {T_0(x),T_1(x),T_2(x)}, where T_0=1, T_1=x, T_2=2x^2-1, i.e.
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x) = b2*x^2 + b1*x + b0,
%    where a = [a2 a1 a0] is sought.
%    Solution:
%      b = [2 0 -2];
%      a = mon2cheb(b);
%
%   See also   MON2CHEB, CHEBPOLY
 
%   Zoltán Csáti
%   2015/03/31
 
%j=1,2 indicates whether we want the resulting polynomial to be written in
%type 1 or 2
 
% Construct the Chebyshev polynomials of the first kind
C = chebpoly(k,j);
% Create the transformation matrix
A = zeros(k+1);
for i = 1:k+1
   A(i:k+1,i) = C{k+2-i};
end

% Perform the basis conversion by solving the linear system
a = A\B(:);
 
end




