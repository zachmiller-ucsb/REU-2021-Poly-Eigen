function varargout = nlevp(name,varargin)
%NLEVP   Collection of nonlinear eigenvalue problems.
%  [COEFFS,FUN,OUT3,OUT4,...] = NLEVP(NAME,ARG1,ARG2,...)
%    generates the matrices and functions defining the problem specified by
%    NAME (a case insensitive string).
%    ARG1, ARG2,... are problem-specific input arguments.
%    All problems are of the form
%      T(lambda)*x = 0
%    where
%      T(lambda)= f0(lambda)*A0 + f1(lambda)*A1 + ... + fk(lambda)*Ak.
%    The matrices A0, A1, ..., Ak are returned in a cell array:
%    COEFFS = {A0,...,Ak}.
%    FUN is a function handle that can be used to evaluate the functions
%    f1(lambda),...,fk(lambda).  For a scalar lambda,
%    F = FUN(lambda) returns a row vector containing
%      F = [f1(lambda), f2(lambda), ..., fk(lambda)].
%    If lambda is a column vector, FUN(lambda) returns a row per element in
%    lambda.
%    [F,FP] = FUN(lambda) also returns the derivatives
%      FP = [f1'(lambda), f2'(lambda), ..., fk'(lambda)].
%    [F,FP,FPP,FPPP,...] = FUN(lambda) also returns higher derivatives.
%    OUT3, OUT4, ... are additional problem-specific output arguments.
%    See the list below for the available problems.
%
%  PROBLEMS = NLEVP('query','problems') or NLEVP QUERY PROBLEMS
%    returns a cell array containing the names of all problems
%    in the collection.
%  NLEVP('help','name') or NLEVP HELP NAME
%    gives additional information on problem NAME, including number and
%    meaning of input/output arguments.
%  NLEVP('query','name') or NLEVP QUERY NAME
%    returns a cell array containing the properties of the problem NAME.
%  PROPERTIES = NLEVP('query','properties') or NLEVP QUERY PROPERTIES
%    returns the properties used to classify problems in the collection.
%  NLEVP('query',property1,property2,...) or NLEVP QUERY PROPERTY1 ...
%    lists the names of all problems having all the specified properties.
%
%  [T,TP,TPP,...] = NLEVP('eval',NAME,LAMBDA,ARG1,ARG2,...)
%  evaluates the matrix function T and its derivatives TP, TPP,...
%  for problem NAME at the scalar LAMBDA.
%
%  NLEVP('version') or NLEVP VERSION
%    prints version, release date, and number of problems
%    of the installed NLEVP collection.
%  V = NLEVP('version')
%    returns a structure V containing version information.
%    V consists of the fields v.number, v.date, and v.problemcount.
%
%  Available problems:
%
%  acoustic_wave_1d   Acoustic wave problem in 1 dimension.
%  acoustic_wave_2d   Acoustic wave problem in 2 dimensions.
%  bicycle            2-by-2 QEP from the Whipple bicycle model.
%  bilby              5-by-5 QEP from Bilby population model.
%  butterfly          Quartic matrix polynomial with T-even structure.
%  cd_player          QEP from model of CD player.
%  closed_loop        2-by-2 QEP associated with closed-loop control system.
%  concrete           Sparse QEP from model of a concrete structure.
%  damped_beam        QEP from simply supported beam damped in the middle.
%  dirac              QEP from Dirac operator.
%  fiber              NEP from fiber optic design.
%  foundation         Sparse QEP from model of machine foundations.
%  gen_hyper2         Hyperbolic QEP constructed from prescribed eigenpairs.
%  gen_tpal2          T-palindromic QEP with prescribed eigenvalues on the
%                     unit circle.
%  gen_tantipal2      T-anti-palindromic QEP with eigenvalues on the unit
%                     circle.
%  gun                NEP from model of a radio-frequency gun cavity.
%  hadeler            NEP due to Hadeler.
%  hospital           QEP from model of Los Angeles Hospital building.
%  intersection       10-by-10 QEP from intersection of three surfaces.
%  loaded_string      REP from finite element model of a loaded vibrating
%                     string.
%  metal_strip        QEP related to stability of electronic model of metal
%                     strip.
%  mirror             Quartic PEP from calibration of cadioptric vision system.
%  mobile_manipulator QEP from model of 2-dimensional 3-link mobile manipulator.
%  omnicam1           9-by-9 QEP from model of omnidirectional camera.
%  omnicam2           15-by-15 QEP from model of omnidirectional camera.
%  orr_sommerfeld     Quartic PEP arising from Orr-Sommerfeld equation.
%  pdde_stability     QEP from stability analysis of discretized PDDE.
%  planar_waveguide   Quartic PEP from planar waveguide.
%  plasma_drift       Cubic PEP arising in Tokamak reactor design.
%  power_plant        8-by-8 QEP from simplified nuclear power plant problem.
%  QEP1               3-by-3 QEP with known eigensystem.
%  QEP2               3-by-3 QEP with known, nontrivial Jordan structure.
%  QEP3               3-by-3 parametrized QEP with known eigensystem.
%  QEP4               3-by-4 QEP with known, nontrivial Jordan structure.
%  QEP5               3-by-3 nonregular QEP with known Smith form.
%  railtrack          QEP from study of vibration of rail tracks.
%  railtrack2         Palindromic QEP from model of rail tracks.
%  relative_pose_5pt  Cubic PEP from relative pose problem in computer vision.
%  relative_pose_6pt  QEP from relative pose problem in computer vision.
%  schrodinger        QEP from Schrodinger operator.
%  shaft              QEP from model of a shaft on bearing supports with a
%                     damper.
%  sign1              QEP from rank-1 perturbation of sign operator.
%  sign2              QEP from rank-1 perturbation of 2*sin(x) + sign(x)
%                     operator.
%  sleeper            QEP modelling a railtrack resting on sleepers.
%  speaker_box        QEP from finite element model of speaker box.
%  spring             QEP from finite element model of damped mass-spring
%                     system.
%  spring_dashpot     QEP from model of spring/dashpot configuration.
%  time_delay         3-by-3 NEP from time-delay system.
%  surveillance       27-by-20 QEP from surveillance camera callibration.
%  wing               3-by-3 QEP from analysis of oscillations of a wing in
%                     an airstream.
%  wiresaw1           Gyroscopic system from vibration analysis of wiresaw.
%  wiresaw2           QEP from vibration analysis of wiresaw with viscous
%                     damping.
%
%  Examples:
%  coeffs = nlevp('railtrack')
%    generates the matrices defining the railtrack problem.
%  nlevp('help','railtrack')
%    prints the help text of the railtrack problem.
%  nlevp('query','railtrack')
%    prints the properties of the railtrack problem.
%
%  For example code to solve all polynomial eigenvalue problems (PEPs)
%  in this collection of dimension < 500 with MATLAB's POLYEIG
%  see NLEVP_EXAMPLE.M.

%  Reference:
%  T. Betcke, N. J. Higham, V. Mehrmann, C. Schroeder, and F. Tisseur.
%  NLEVP: A Collection of Nonlinear Eigenvalue Problems,
%  MIMS EPrint 2011.116, Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, 2011

% Check inputs
if nargin < 1, error('Not enough input arguments'); end
if ~ischar(name), error('NAME must be a string'); end

name = lower(name);

if strcmp(name,'query')
    if nargin == 1
       error('Not enough input arguments')
    end
    [varargout{1:nargout}] = nlevp_query(varargin{:});
    return
end

if strcmp('string',name)
    name = 'spring';
    warning('NLEVP:string_renamed','Problem string has been renamed spring.')
end

if strcmp('version',name)
    [varargout{1:nargout}] = nlevp_version(varargin{:});
    return
end

switch name
    case 'help'
        if nargin < 2
           help nlevp
        else
           if ~nlevp_isoctave
               help(varargin{1})
           else
               % Uglier code necessary for Octave.
               eval(['help ', varargin{1}]);
           end
        end
    case 'eval'
        [varargout{1:max(nargout,1)}] = nlevp_eval(varargin{:});
    otherwise
        [varargout{1:nargout}] = feval(name,varargin{:});
end
