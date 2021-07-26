function matrix = evgen_normal1_center(n, ~, ~)
% Return eigenvalues taken from a normal distribution with mean zero and
% variance 0.2, clamped to three standard deviations.

    matrix = stdrandncutoff(1, n, 0.6);
    matrix = sort(matrix);
end