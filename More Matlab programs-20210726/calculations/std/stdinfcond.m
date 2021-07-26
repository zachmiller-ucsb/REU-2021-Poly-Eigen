function answer = stdinfcond(M)
% Return the infinity condition number of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = cond(numerizedM, Inf);
end
