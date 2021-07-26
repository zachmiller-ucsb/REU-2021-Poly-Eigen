function answer = stdtwonorm(M)
% Return the two norm of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = norm(numerizedM);
end
