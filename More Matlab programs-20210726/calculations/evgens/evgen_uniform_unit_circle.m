function matrix = evgen_uniform_unit_circle(n, ~, ~)
% Return n eigenvalues from a uniform distribution in the node range.  Note
% that if the node range is a purely real or imaginary line, the node range
% is interpreted to be a line.  Otherwise, the node range is interpreted to
% be a disk in the complex plane.

    noderange = [-1-1i 1+1i] / sqrt(2);
    matrix = stdrandnmtx(1, n, noderange);
end
