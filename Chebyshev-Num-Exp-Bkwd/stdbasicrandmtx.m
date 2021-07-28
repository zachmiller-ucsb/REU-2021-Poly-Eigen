function mtx = stdbasicrandmtx(r, c, range, stdrandfn)
% Return a matrix with scalar coefficients from the given distribution
% within the range.  If the range is a purely real or imaginary line, the
% range is interpreted to be a line.  Otherwise, the node range is
% interpreted to be a disk in the complex plane.  Use the given stdrandfn
% random function.
%
% Do not call this function directly, call stdrandnmtx() or stdrandumtx()
% instead.

    [~, ~, mean, radius] = stdexplainrange(range);
    if imag(radius) == 0
        mtx = stdrandfn(r, c, real(radius)) + mean;
        return;
    end
    if real(radius) == 0
        mtx = stdrandfn(r, c, imag(radius)) + mean;
        mtx = mtx * 1i;
        return;
    end
    radii = sqrt(stdrandfn(r, c, 1));
    thetas = stdrandfn(r, c, 2*pi);
    thetacosines = cos(thetas);
    thetasines = sin(thetas);
    mtx = thetacosines .* radii + thetasines .* radii * 1i;
    mtx = mtx * radius + mean;
end
