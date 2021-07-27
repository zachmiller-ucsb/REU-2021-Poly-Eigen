function [coeffs,fun] = speaker_box
%SPEAKER_BOX  QEP from model of a speaker box.
%  [COEFFS,FUN] = nlevP('speaker_box') generates the coefficient matrices of
%  a quadratic matrix polynomial lambda^2*A2 + lambda*A1 + A0 arising from
%  a model of a speaker box.
%  This problem is real symmetric and has dimension 107.
%  The matrices are sparse and returned in a cell array:
%  COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, real, symmetric.

load('speaker_box.mat');

coeffs = {C, B, A};
fun = @(lam) nlevp_monomials(lam,2);
