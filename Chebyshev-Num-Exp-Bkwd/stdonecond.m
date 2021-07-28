function answer = stdonecond(M)
% Return the 1-condition number of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = cond(numerizedM, 1);
end
