function [kratios, bounds] = ev_condition_chebyshev(d, polysize, ...
    nodegen, polygen, evgen, pkmean, pkwidth, boundgens)
% Chebyshev eigenvalue condition number bound test function.
% Use this function from the generic_ev_test script.

    if ~iscell(nodegen) || ~isequal(nodegen{1}, @nodes_chebyshev_d_f)
        error 'This experiment can only run with the Chebyshev nodes.'
    end
    [kratios, bounds] = ev_generic_condition( ...
        @std_xchg_chebyshev_monomial, @basiseval_chebyshev, ...
        d, polysize, nodegen, polygen, evgen, pkmean, pkwidth, boundgens);
end
