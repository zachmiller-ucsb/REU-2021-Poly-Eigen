function matrix = evgen_normal1_left_tail(n, ~, ~)
% Return eigenvalues taken from a normal distribution with mean -1 and
% variance 0.2, clamped to three standard deviations, and no less than -1.

    matrix = abs(evgen_normal1_center(n, [], [])) - 1;
    matrix = sort(matrix);
end