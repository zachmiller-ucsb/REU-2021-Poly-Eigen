function [rangemin, rangemax, center, radius] = stdexplainrange(noderange)
% Return the values above describing the general shape of noderange.

    rangemin = min(real(noderange)) + min(imag(noderange)) * 1i;
    rangemax = max(real(noderange)) + max(imag(noderange)) * 1i;
    center = (rangemax + rangemin) / 2;
    radius = (rangemax - rangemin) / 2;
end
