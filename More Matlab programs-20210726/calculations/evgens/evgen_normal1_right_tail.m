function matrix = evgen_normal1_right_tail(n, ~, ~)
% Return eigenvalues taken from a normal distribution with mean 1 and
% variance 0.2, clamped to three standard deviations, and no more than 1.

    matrix = 1 - abs(evgen_normal1_center(n, [], []));
    matrix = sort(matrix, 'descend');
end