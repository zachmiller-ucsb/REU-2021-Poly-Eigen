function hash = stdrandomhash(d, polysize, pkmean, pkwidth)
% Return a deterministic hash of the given inputs.

    hash = d;
    hash = hashwith(hash, polysize);
    hash = hashwith(hash, pkmean);
    hash = hashwith(hash, pkwidth);
    hash = uint32(hash);
end


function hash = hashwith(state, input)
% A multiplicative hash function with factor 1664525.  This is reasonable
% here because the purpose is to produce a random number generator seed.

    hash = mod(state * 1664525 + input, 2^32);
end
