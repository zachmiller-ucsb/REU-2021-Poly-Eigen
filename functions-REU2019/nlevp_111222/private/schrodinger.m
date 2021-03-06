function [coeffs,fun] = schrodinger
%SCHRODINGER  QEP from Schrodinger operator.
%  [COEFFS,FUN] = nlevp('schrodinger') generates the coefficient matrices of
%  a quadratic matrix polynomial lambda^2*A + lambda*B + C of dimension
%  1998.  The spectrum of this matrix polynomial is the second order
%  spectrum of the Schrodinger operator with nonzero potential,
%  relative to the subspace generated by fourth order Hermite elements
%  on a uniform mesh on the interval [-49,49].
%  The matrices are returned in a cell array: COEFFS = {C, B, A}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, real, symmetric, sparse.

%  Reference:
%  L. Boulton and M. Levitin, On approximation of the eigenvalues of
%  perturbed periodic Schroedinger operators, J. Phys. A: Math. Theor. 40
%  (2007), pp. 9319-9329

load schrodinger
coeffs = C;
% Symmetrize (already symmetric up to roundoff).
for i = 1:3
    coeffs{i} = (coeffs{i} + coeffs{i}')/2;
end
fun = @(lam) nlevp_monomials(lam,2);
