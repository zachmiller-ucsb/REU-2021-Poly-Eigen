function answer = stdcond(M)
% Return the two norm condition number of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = cond(numerizedM);
end
