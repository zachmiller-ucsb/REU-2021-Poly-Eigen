function [coeffs,fun] = railtrack
%RAILTRACK   QEP from study of vibration of rail tracks.
%  [COEFFS,FUN] = nlevp('railtrack') generates the coefficient matrices
%  of the railtrack problem, a quadratic eigenvalue problem
%      (lambda^2 C + lambda B + A)*x = 0,
%  which models vibration of rail tracks under the excitation of high speed
%  trains.  This problem is T-palindromic, i.e., C = A.' and B = B.',
%  and has dimension 1005.
%  The matrices are sparse and returned in a cell array: COEFFS = {A, B, C}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, t-palindromic, sparse.

% References:
% A. Hilliges. Numerische L\"osung von quadratischen Eigenwertproblemen
%    mit Anwendung in der Schienendynamik. Master's thesis, TU Berlin, 2004.
% Andreas Hilliges, Christian Mehl, and Volker Mehrmann. On the
%    solution of palindromic eigenvalue problems. In P. Neittaanmaki,
%    T. Rossi, S. Korotov, E. Onate, J. Periaux, and D. Knorzer,
%    editors, Proceedings of the European Congress on Computational
%    Methods in Applied Sciences and Engineering (ECCOMAS 2004),
%    Jyvaskyla, Finland, 2004.
%    http://www.mit.jyu.fi/eccomas2004/proceedings/proceed.html
% D. Steven Mackey, Niloufer Mackey, Christian Mehl, and Volker Mehrmann.
%    Structured polynomial eigenvalue problems: Good vibrations from
%    good linearizations. SIAM J. Matrix Anal. Appl., 28(4):1029-1051, 2006.

load('railtrack_A.mat') % loads sA
load('railtrack_B.mat') % loads sB

coeffs = {sA.', sB, sA};

fun = @(lam) nlevp_monomials(lam,2);
