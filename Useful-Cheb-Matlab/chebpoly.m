function C = chebpoly(k,j)
% C = CHEBPOLY(N) returns Chebyshev polynomials of the first or second kind from 
% degree 0 to n
C = cell(k+1,1);   % using cell array because of the different matrix sizes
C{1} = sym(1, 'f');  % T_0 = 1
C{2} = [sym(j,'f'), sym( 0, 'f')]; % T_1 = x ; T_2 = 2x
 
for i = 3:k+1      % T_n = 2xT_n-1 - T_n-2
    C{i} = [2*C{i-1} 0] - [0 0 C{i-2}];
end
end
