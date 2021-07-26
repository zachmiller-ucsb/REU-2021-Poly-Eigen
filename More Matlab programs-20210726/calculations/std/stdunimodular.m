function mtx = stdunimodular(r, c, mean, radius)
% Return a random unimodular matrix, with coefficients scaled to the given
% radius around the provided mean.

    mtxcond = stdcondmax();
    mtx = zeros(r, c);
    while mtxcond >= stdcondmax()
        mtx = stdrandnmtx(r, c, [mean-radius, mean+radius]);
        mtxcond = cond(mtx);  % not stdcond(), the precision is known here
    end
end
