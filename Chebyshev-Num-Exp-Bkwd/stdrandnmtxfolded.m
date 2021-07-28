function mtx = stdrandnmtxfolded(r, c, range)
% Return a matrix with scalar coefficients from a normal distribution
% within the range.  If the range is a purely real or imaginary line, the
% range is interpreted to be a line.  Otherwise, the node range is
% interpreted to be a disk in the complex plane.

    mtx = stdbasicrandmtx(r, c, range, @stdrandnfolded);
end
