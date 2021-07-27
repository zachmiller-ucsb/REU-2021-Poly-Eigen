function [coeffs,fun] = shaft
%SHAFT  QEP from model of a shaft on bearing supports with a damper.
%  [COEFFS,FUN] = nlevp('shaft') generates the coefficient matrices of a
%  quadratic matrix polynomial lambda^2*A2 + lambda*A1 + A0 arising from
%  a finite element model of a shaft on bearing supports with a damper.
%  This problem is real symmetric and has dimension 400.
%  The matrices are sparse and returned in a cell array:
%  COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, real, symmetric, sparse.

load('shaft.mat');
coeffs = {K, C, M};
fun = @(lam)nlevp_monomials(lam,2);
