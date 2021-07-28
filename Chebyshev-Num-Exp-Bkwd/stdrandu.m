function sample = stdrandu(r, c, radius)
% Return a sample from the normal distribution, scaled to the given radius.

    sample = unifrnd(-1, 1, [r c]) * radius;
end
