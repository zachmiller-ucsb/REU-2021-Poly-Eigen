function [kratios, bounds] = ev_generic_condition(xchggen, basiseval, ...
    d, polysize, nodegen, polygen, evgen, pkmean, pkwidth, boundgens)
% Generic eigenvalue condition number bound test function.
% Do not use this function from the generic_ev_test script.

    preamble_numeric();

    rndstate = stdstartdeterministicrandom(d, polysize, pkmean, pkwidth);

    % obtain nodes, eigenvalues, and matrix polynomial coefficients Pk
    [nodes, noderange] = nodegeneval(nodegen, d);
    [evs, pks] = polygen(d, polysize, ...
        nodes, noderange, evgen, pkmean, pkwidth);

    % obtain corresponding matrix polynomial coefficients Bk
    x2m = xchggen(d, nodes, noderange);
    bks = polygen_rewrite_polynomial(x2m, pks);

    % evaluate the eigenvalue ratios
    kratios = kratioseval(d, nodes, noderange, ...
        evs, pks, bks, basiseval);

    % evaluate the bounds
    bounds = boundseval(boundgens, d, nodes, noderange, evs, x2m);

    stdstopdeterministicrandom(rndstate);
end
