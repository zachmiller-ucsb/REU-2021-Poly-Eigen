function answer = stdtwonormmax(M)
% Return an upper bound on the two norm of M.  Use that the two norm of M
% is bounded above by the sqrt of the product of M's one and inf norms,
% and that the one norm of an nxn matrix is less than n times its inf norm.
% In total, the two norm is bounded above by sqrt(n) * infnorm(M).  Allow
% for a (rather very pessimistic) 0.1% relative error due to implicit
% numerization in stdinfnorm().

    answer = sqrt(length(M)) * stdinfnorm(M) * 1.001;
end
