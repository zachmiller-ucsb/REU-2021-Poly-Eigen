function sample = stdrandnfolded(r, c, radius)
% Return a sample from the normal distribution, scaled to the given radius.
% Fold the tails beyond stdnormalmax() standard deviations.

    normalmax = stdnormalmax();
    sample = randn(r, c);
    sample = max(-normalmax, sample);
    sample = min(sample, normalmax);
    sample = sample * radius / normalmax;
end
