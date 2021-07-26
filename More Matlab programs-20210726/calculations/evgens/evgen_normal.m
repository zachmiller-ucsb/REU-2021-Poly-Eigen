function matrix = evgen_normal(n, ~, noderange)
% Return n eigenvalues from a normal distribution in the node range.  Note
% that if the node range is a purely real or imaginary line, the node range
% is interpreted to be a line.  Otherwise, the node range is interpreted to
% be a disk in the complex plane.

    matrix = stdrandnmtx(1, n, noderange);
end
