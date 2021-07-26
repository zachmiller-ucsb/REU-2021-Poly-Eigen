function sample = stdrandncutoff(r, c, radius)
% Return a sample from the normal distribution, scaled to the given radius.
% Cut off the tails beyond stdnormalmax() standard deviations.

    normalmax = stdnormalmax();
    sample = randn(r, c);
    for k = 1:r
        for j = 1:c
            while abs(sample(k, j)) > normalmax
                sample(k, j) = normrnd(0, 1);
            end
        end
    end
    sample = sample * radius / normalmax;
end
