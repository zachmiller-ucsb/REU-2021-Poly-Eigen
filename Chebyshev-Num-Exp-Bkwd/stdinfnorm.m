function answer = stdinfnorm(M)
% Return the infinity norm of the (possibly symbolic) matrix M.

    numerizedM = stdnumerize(M);
    answer = norm(numerizedM, Inf);
end
