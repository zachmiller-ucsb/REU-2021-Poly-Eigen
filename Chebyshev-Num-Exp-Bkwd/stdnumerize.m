function answer = stdnumerize(sym)
% For the purposes of this program, computing matrix norms and condition
% numbers can be done in VPA arithmetic, rather than symbolically.  This
% function is used to do that, and also clearly marks where that conversion
% happens.

    answer = vpa(sym);
end
