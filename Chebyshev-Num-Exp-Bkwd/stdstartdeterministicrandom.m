function rndstate = stdstartdeterministicrandom(d, polysize, pkmean, pkwidth)
% Set the random number generator for deterministic operation, with a seed
% constructed from the given arguments.

    rndstate = rng;
    newseed = stdrandomhash(d, polysize, pkmean, pkwidth);
    % must specify the rng type because the local and parallel default rng
    % types are different for some reason :(.
    rng(newseed, 'twister');
end
