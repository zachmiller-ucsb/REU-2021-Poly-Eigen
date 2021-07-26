function matrix = evgen_uniform(n, ~, noderange)
% Return n eigenvalues from a uniform distribution in the node range.  Note
% that if the node range is a purely real or imaginary line, the node range
% is interpreted to be a line.  Otherwise, the node range is interpreted to
% be a disk in the complex plane.

    matrix = stdrandumtx(1, n, noderange);
end
