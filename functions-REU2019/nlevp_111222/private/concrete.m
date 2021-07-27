function [coeffs,fun] = concrete(mu)
%CONCRETE  Sparse QEP from model of a concrete structure.
%  [COEFFS,FUN] = nlevp('concrete',MU) generates the coefficient matrices
%  of a quadratic matrix polynomial lambda^2*C + lambda*B + A arising in a
%  model of a concrete structure supporting a machine assembly.
%  This problem is complex symmetric and has dimension 2472.
%  The matrices are sparse and returned in a
%  cell array: COEFFS = {A, B, C}.
%  C = mass matrix, real diagonal, low rank.
%  B = viscous damping matrix, purely imaginary
%      and diagonal, B = i*C_v, low rank.
%  A = stiffness + uniform hysteretic damping
%      matrix, A = (1 + i*MU)*K.
%  By default, MU = 0.04.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, symmetric, parameter-dependent,
%  sparse.

%  Reference:
%  A. Feriani, F. Perotti and V. Simoncini,
%  Iterative system solvers for the frequency analysis of
%  linear mechanical systems,
%  Computer Methods in Appl. Mech. Eng., 190 (2000), pp. 1719-1739.

if nargin < 1, mu = 0.04; end;

load('concrete_K.mat') % loads K_c1
load('concrete_C.mat') % loads C_c1
load('concrete_M.mat') % loads M_c1

coeffs = {(1+mu*1i)*K_c1, C_c1, M_c1};

fun = @(lam) nlevp_monomials(lam,2);
