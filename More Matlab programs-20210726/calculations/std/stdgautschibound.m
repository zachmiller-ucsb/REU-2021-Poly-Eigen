function answer = stdgautschibound(d, nodes)
% Return the Gautschi bound on the infinity norm of the inverse of the
% Vandermonde matrix associated with the given nodes.

    g = ones(1, d+1);
    for k = 1:d+1
        for j = [1:k-1 k+1:d+1]
            g(k) = g(k) * ((1+abs(nodes(j))) / abs(nodes(k)-nodes(j)));
        end
    end
    answer = max(g);
end
