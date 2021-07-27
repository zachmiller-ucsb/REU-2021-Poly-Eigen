function [coeffs,fun] = hospital
%HOSPITAL   QEP from model of Los Angeles Hospital building.
%  [COEFFS,FUN] = nlevp('hospital') constructs a 24-by-24 quadratic
%  matrix polynomial lambda^2*M + lambda*D + K arising in a model of the
%  Los Angeles University Hospital building.  There are 8 floors,
%  each with 3 degrees of freedom.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, real.

%  Reference:
%  Y. Chahlaoui and P. M. Van Dooren. A collection of benchmark examples
%  for model reduction of linear time invariant dynamical systems.
%  MIMS EPrint 2008.22, Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, 2002.

load hospital_KD.mat

n = size(D,1);
M = eye(n);

coeffs = {K, D, M};

fun = @(lam) nlevp_monomials(lam,2);
