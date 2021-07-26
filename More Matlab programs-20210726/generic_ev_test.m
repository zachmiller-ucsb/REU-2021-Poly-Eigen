%%%%%%%%%%%%%%%%%%%% Generic test script.

% Use this file as a template to run experiments.  However, do not modify
% this file, so that the template stays clean at all times (and so that
% particular experiments can be version controlled independently).

% This script uses the Parallel Toolbox to speed up computation.  Make sure
% it is installed, and that you have a local parallel pool ready to go.

% Run the preamble, which sets the running conditions for the scripts.
preamble;

% Set the path to save figures.  No trailing slashes.  Defaults to the
% current MATLAB folder.
figuresavepath = '.';

% Configure the figures.
plotfontsize = 12;
plottitle = true;

%%%%%%%%%%%%%%%%%%%% Set up conditions.

% Matrix polynomial size.  Here, matrix polynomials are expected to be
% regular, hence they must be square: nxn.
polysize = 2;

% Matrix polynomial degrees.  The bounds on eigenvalue conditioning depend
% strongly on the degree, d.  So, provide a range of degrees to run
% experiments with.
degreerange = 1:12;

% Choice of interpolation basis.  This choice affects other parameters, see
% those in the corresponding files.  Pick one.

experiment = @ev_condition_chebyshev;
% experiment = @ev_condition_lagrange;
% experiment = @ev_condition_newton;
% experiment = @ev_condition_taylor;


% Choice of bounds.  Select the bounds to plot corresponding to the
% experiment.
boundgens = {
%    @bounds_generic
% Chebyshev specific
%    @bounds_chebyshev_dp1_evs
%    @bounds_chebyshev_evs
%    @bounds_chebyshev_norm
% Chebyshev specific --- unused
%    @bounds_chebyshev_cond_2
%    @bounds_chebyshev_upperbounds
% Lagrange specific
%    @bounds_lagrange_g
%    @bounds_lagrange_l
% Lagrange specific --- only valid for roots of unity
%    @bounds_lagrange_lru
% Lagrange specific --- only valid for roots of unity on unit circle
%    @bounds_lagrange_lru_unit_circle
% Lagrange specific --- only valid for Leja sequences
%    @bounds_lagrange_leja
% Lagrange specific --- unused
%    @bounds_lagrange
%    @bounds_lagrange_d
% Lagrange specific --- best choice of d and g variants above
%    @bounds_lagrange_dg
% Newton specific
%    @bounds_newton
};


% Choice of interpolation node interval / region.  Some examples are listed
% below for convenience.  For Chebyshev experiments, the nodes are expected
% to be real, and between zero and one.  These node functions will be
% called by the experiment for each degree.  Select one.

% Chebyshev nodes.  For these generators, provide a scaling factor as a
% parameter.  Note that the Chebyshev experiment can only use the real
% Chebyshev nodes.
nodegen = {@nodes_chebyshev_d_f, [1]};
% nodegen = {@nodes_chebyshev_d_f, [1i]};

% Complex circles based on the roots of unity.  For these generators,
% provide the center and the radius as parameters.
% nodegen = {@nodes_unity_roots_d_cr, [-5e-1, 1e-4]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1e-2]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [0, 1e4]};
% nodegen = {@nodes_unity_roots_d_cr, [5e-1, 1e-2]};
% nodegen = {@nodes_unity_roots_d_cr, [5e-1, 5e-1]};
% nodegen = {@nodes_unity_roots_d_cr, [1, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1/2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [1e1, 1e4]};
% nodegen = {@nodes_unity_roots_d_cr, [1e2, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e2, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e2, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e2, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [1e2, 1e4]};
% nodegen = {@nodes_unity_roots_d_cr, [1e3, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e3, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e3, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e3, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [1e3, 1e4]};
% nodegen = {@nodes_unity_roots_d_cr, [1e4, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e4, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e4, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e4, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [1e4, 1e4]};
% nodegen = {@nodes_unity_roots_d_cr, [1e5, 1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e5, 1e1]};
% nodegen = {@nodes_unity_roots_d_cr, [1e5, 1e2]};
% nodegen = {@nodes_unity_roots_d_cr, [1e5, 1e3]};
% nodegen = {@nodes_unity_roots_d_cr, [1e5, 1e4]};

% Complex circles based on the roots of unity.  For these generators,
% provide the center as a parameter (the radius is 1 / (d+1)).
% nodegen = {@nodes_unity_roots_d_c_1odp1, [0]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1e1]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1e2]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1e3]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1e4]};
% nodegen = {@nodes_unity_roots_d_c_1odp1, [1e5]};

% Complex circles based on the roots of unity.  For these generators,
% provide the center as a parameter (the radius is based on the Leja upper
% bound).
% nodegen = {@nodes_unity_roots_d_c_lejaupper, [1e1]);

% Complex circles based on the roots of unity.  For these generators,
% provide the center as a parameter (the radius is sin(pi/(d+1))).
% nodegen = {@nodes_unity_roots_d_c_sinpiodp1, [0]);

% Complex disks.  For these generators, provide the center and the radius
% as parameters.  These nodes are normally distributed.
% nodegen = {@nodes_disk_d_cr, [0 1/4]};
% nodegen = {@nodes_disk_d_cr, [0 1]};
% nodegen = {@nodes_disk_d_cr, [0 1e2]};
% nodegen = {@nodes_disk_d_cr, [1e1 1/2]};
% nodegen = {@nodes_disk_d_cr, [1e1 5]};
% nodegen = {@nodes_disk_d_cr, [1e2 1/4]};
% nodegen = {@nodes_disk_d_cr, [1e2 1]};
% nodegen = {@nodes_disk_d_cr, [1e2 1e2]};
% nodegen = {@nodes_disk_d_cr, [1e5 1]};
% nodegen = {@nodes_disk_d_cr, [1e5 1e2])};

% Complex disks.  For these generators, provide the center and the radius
% as parameters.  These nodes are uniformally distributed.
% nodegen = {@nodes_udisk_d_cr, [0 1]};
% nodegen = {@nodes_udisk_d_cr, [1e2 1e2]};

% Equidistant nodes.  For these generators, provide the start node, the end
% node, and a factor to multiply them all as parameters.
% nodegen = {@nodes_equidistant_d_abf, [-1 0 1]};
% nodegen = {@nodes_equidistant_d_abf, [-1 0 1i]};
% nodegen = {@nodes_equidistant_d_abf, [0 1 1]};
% nodegen = {@nodes_equidistant_d_abf, [0 1 1i]};
% nodegen = {@nodes_equidistant_d_abf, [-1 1 1]};
% nodegen = {@nodes_equidistant_d_abf, [-1 1 1i]};

% Equidistant nodes between 1 and d+1, times a scaling factor.  For these
% generators, provide the factor as a parameter.
% nodegen = {@nodes_equidistant_d_1dp1_f, [1]};
% nodegen = {@nodes_equidistant_d_1dp1_f, [1i]};

% Equidistant nodes between -d and d, times a scaling factor.  For these
% generators, provide the factor as a parameter.
% nodegen = {@nodes_equidistant_d_pmd_f, [1]};
% nodegen = {@nodes_equidistant_d_pmd_f, [1i]};

% Exponential nodes with assorted bases.  For these generators, pass the
% base of the exponential and an overall factor as parameters.
% nodegen = {@nodes_expplus_d_bf, [(11/10) 1]};
% nodegen = {@nodes_expplus_d_bf, [(3/2) 1]};
% nodegen = {@nodes_expplus_d_bf, [3 1]};

% nodegen = {@nodes_expplus_d_bf, [2 1]};
% nodegen = {@nodes_expplus_d_bf, [2 1i]};
% nodegen = {@nodes_expplus_d_bf, [2 -1]};
% nodegen = {@nodes_expplus_d_bf, [2 -1i]};

% nodegen = {@nodes_expplus_d_bf, [(1+1i) 1]};
% nodegen = {@nodes_expplus_d_bf, [(2+2i) 1]};

% nodegen = {@nodes_expminus_d_bf, [2 1]};
% nodegen = {@nodes_expminus_d_bf, [2 1i]};
% nodegen = {@nodes_expminus_d_bf, [2 -1]};
% nodegen = {@nodes_expminus_d_bf, [2 -1i]};

% Exponential nodes, two sided.
% nodegen = {@nodes_expplus_twosided_d_bf, [2 1]};
% nodegen = {@nodes_expminus_twosided_d_bf, [2 1]};
% nodegen = {@nodes_expplusminus_twosided_d_bf, [2 1]};


% Choice of polynomial generator.  Pick one.

% This generator is limited because the polynomial matrix coefficients
% will be singular except for the first one.
% polygen = @polygen_smith;

% This generator (virtually always) produces matrix coefficients that are
% regular scalar matrices.  This is the default.
polygen = @polygen_pseudosmith;


% Choice of eigenvalue generator.  Pick one.

% Eigenvalues distributed in the noderange with uniform distribution.  This
% is the default.
evgen = @evgen_uniform;

% Eigenvalues distributed in the noderange with normal distribution.
% evgen = @evgen_normal;

% Eigenvalues with normal distribution between -1 and 1, and selected
% accumulation points.
% evgen = @evgen_normal1_center;
% evgen = @evgen_normal1_left_tail;
% evgen = @evgen_normal1_right_tail;
% evgen = @evgen_normal1_two_tail;

% Eigenvalues distributed between the origin and the node with smallest
% absolute value, bunched closer to said node in proportion to the
% square root of the linear fraction.
% evgen = @evgen_evtoxl90pc;

% Eigenvalues distributed towards the center of the node range.
% evgen = @evgen_evtocenter90pc;

% Eigenvalues distributed equidistantly between the center of the
% node range and the node with largest absolute value.
% evgen = @evgen_evtoeq;

% Eigenvalues distributed between the center of the node range and the node
% with largest largest absolute value, bunched closer to said node in
% proportion to the square root of the linear fraction.
% evgen = @evgen_evtoxh;
% evgen = @evgen_evtoxh90pc;

% Eigenvalues distributed near the node with largest absolute value, and
% the node furthest away from that node.
% evgen = @evgen_evdipole90pc;
% evgen = @evgen_evdipole99pc;

% Eigenvalues for fixed areas --- use judiciously.
% evgen = @evgen_uniform_unit_circle;

% Eigenvalues being the roots of unity, distributed on a circle with radius
% set to the smallest nonzero absolute value of the nodes.  The second
% generator scales the radius even lower by a factor of 10^-d.
% evgen = @evgen_unity_roots_0_xl;
% evgen = @evgen_unity_roots_0_xlemn;


% Choice of mean and width of normal distribution controlling the
% polynomial matrix coefficient generation.
pkmean = 0;
pkwidth = 50;


%%%%%%%%%%%%%%%%%%%% Set up result cell arrays.

% Meaning of each dimension of bound_data:
%   first level is per bound generator
%   second level is per upper / lower bound
%   third level is per degree
bound_data = zeros(length(boundgens), 2, max(degreerange));

% Meaning: scatter data for the eigenvalue condition number ratios.
ratio_data = cell(max(degreerange), 1);

%%%%%%%%%%%%%%%%%%%% Parallel for loop for actual computation.

% Due to parfor index restrictions, will have to swizzle results later.
par_ratio_data = cell(max(degreerange), 1);
par_bound_data = cell(max(degreerange), 1);

% Reverse the degree indexing so that parallel processing schedules better.
% Unfortunately, parfor requires a strictly increasing integer range.
dmax = length(degreerange);

% Known complaints from the MATLAB linter.
%   * That degree is a broadcast variable --- negligible expected impact.

tic
parfor dk = 1:dmax
    d = degreerange(dmax - dk + 1);
    [kratios, bounds] = experiment( ...
        d, polysize, nodegen, polygen, evgen, pkmean, pkwidth, boundgens);
    par_ratio_data{dk} = kratios;
    par_bound_data{dk} = bounds;
end
toc

% Swizzle results back into shape.
for dk = 1:dmax
    dwrite = degreerange(dmax - dk + 1);
    ratio_data(dwrite) = par_ratio_data(dk);
    for j = 1:length(boundgens)
        bound_data(j, 1, dwrite) = par_bound_data{dk}(j, 1);
        bound_data(j, 2, dwrite) = par_bound_data{dk}(j, 2);
    end
end

%%%%%%%%%%%%%%%%%%%% Plotting of the data.

% Plotting method to display results.  Pick one.

% This plotting method is intended for development and interactive work.
% It opens a figure with all bounds displayed.
plotgen = @generic_plot;

% This plotting method is intended for batch operation.  It will
% automatically save the plots as png and fig files, and do so with all
% possible combinations of provided bounds.
% plotgen = @generic_autosave_plot;

% Eigenvalue scatterplot colorizers.  Pick one.
plotscattercolorizer = @plotscatter_black;
% plotscattercolorizer = @plotscatter_grad_blue;
% plotscattercolorizer = @plotscatter_grad_blue2x;
% plotscattercolorizer = @plotscatter_grad_redgreen;

% Actually plot the results.
plotgen(plotscattercolorizer, plotfontsize, plottitle, ...
    experiment, degreerange, nodegen, polygen, evgen, ...
    boundgens, bound_data, ratio_data, figuresavepath);
