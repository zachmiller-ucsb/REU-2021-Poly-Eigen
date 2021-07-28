function answer = stdonenorm(M)
% Return the one norm of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = norm(numerizedM, 1);
end
