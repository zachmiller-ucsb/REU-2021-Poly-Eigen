function V = stdvander(nodes)
% Return the standard form of the Vandermonde matrix we use here.

    V = fliplr(vander(nodes));
end
