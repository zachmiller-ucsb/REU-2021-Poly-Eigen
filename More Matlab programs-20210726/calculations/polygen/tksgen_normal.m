function answer = tksgen_normal(d, polysize, noderange)
% Return values for the function T, assuming T is polysize x polysize, and
% that there will be d+1 interpolation nodes.

    answer = stdrandnmtx((d+1) * polysize, polysize, noderange);
end
