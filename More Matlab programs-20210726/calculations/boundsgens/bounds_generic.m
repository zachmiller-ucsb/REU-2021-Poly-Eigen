function [u, l] = bounds_generic(~, ~, ~, ~, L)
% Calculate the generic bounds for the ratio of condition numbers.

    u = stdonecond(L);
    l = 1/u;
end
