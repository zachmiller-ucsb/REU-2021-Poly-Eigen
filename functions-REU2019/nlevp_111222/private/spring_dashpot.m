function [coeffs,fun] = spring_dashpot(n,rho,state)
%SPRING_DASHPOT  QEP from model of spring/dashpot configuration.
%  [COEFFS,FUN] = nlevp('spring_dashpot',N,RHO,S) generates an
%  N-by-N quadratic matrix polynomial
%  Q(lambda) = lambda^2*M + lambda*D + K
%  arising from the modelling of a linear spring in parallel with
%  N/2-1 Maxwell elements, where a Maxwell element is a spring in series
%  with a dashpot.  RHO is the material density.
%  N should be even; otherwise N+1 is used.
%  This example reflects the structure only, since the matrices themselves
%  are not from a finite element model but randomly generated to have
%  the desired properties of symmetry etc.
%  S (default: S = 0) is used in the function call RAND('twister',S) to
%  set the random number generator (see HELP RAND); if S < 0 then
%  the newer syntax RNG(-S) is used.
%  Set S = 'noseed' to leave the random number generator unseeded.
%  The singular mass matrix M contributes infinite eigenvalues.
%  By default, N = 10 and RHO = 7850.
%  Properties of the coefficient matrices:
%     M is rank deficient and symmetric.
%     D is rank deficient and block diagonal.
%     K is symmetric and has arrowhead structure.
%  The matrices are returned in a cell array: coeffs = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties qep, pep, real, parameter-dependent,
%  scalable, sparse, random.

%  Reference:
%  Anthony Gotts, Report regarding model reduction, model compaction
%  project, University of Nottingham, Feb. 2005.

%  Default is steel for rho, the material density.
if nargin < 1 || isempty(n)
    n = 10;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
           'now targets the true dimension of the problem']);
end
if nargin >= 1
    warning('NLEVP:argument_change',['Note that the first input ',...
           'argument is now the problem dimension.']);
end
if nargin < 2 || isempty(rho), rho = 7850; end
if nargin < 3 || isempty(state), state = 0; end

if ~strcmpi(state,'noseed')
    if state >= 0
        rand('twister',state) % Set state for random numbers.
    else
        rng(-state)           % New syntax.
    end
end

m = round(n/2)-1;
n = 2*(m+1);

% Generate one random element stiffness matrix and reuse.
K = rand(2); K = K*K';

% Generate one random element mass matrix and reuse.
M = rand(2); M = M*M';
e = rand(1,m);
eta = rand(1,m);
alpha_rho = sum(e);

B = kron(eta, K);

% M = blkdiag(rho*M, zeros(2*m));
M = sparse([1;2;1;2],[1;1;2;2],rho*M(:),n,n);

% D = kron(blkdiag(0,diag(eta)),K);
D = kron(diag(sparse([0,eta])),K);

% K = kron(blkdiag(alpha_rho,diag(e)), K) + [zeros(2) B; B' zeros(2*m)];
K = kron(diag(sparse([alpha_rho,e])), K);
K(1:2,3:end) = B;
K(3:end,1:2) = B';

coeffs{3} = M;
coeffs{2} = D;
coeffs{1} = K;
fun = @(lam) nlevp_monomials(lam,2);
