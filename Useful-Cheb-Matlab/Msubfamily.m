function [M1,M0] = Msubfamily(coeff,ep,j)

% Constructs body of subfamily linearization for a matrix polynomial P
% written in the chebyshev basis dependent on parameters 0 <= ep <= k - 1 (size), k
% (degree of P), and j for chebyshev basis (type 1 --> j = 1, type 2 --> j = 2)

% Matrix coefficients A_i of P stored in block vector coeff

L = length(coeff); % Matrices are 1-indexed so matrix polynomial will be of degree L - 1, i.e. A_i = coeff(i + 1)
k = L - 1;
nTemp = coeff(1); % Temp used to obtain size of P
n = nTemp(1); % Size of P
M1 = sym(zeros((k-ep)*n, (ep+1)*n));
M0 = sym(zeros((k-ep)*n, (ep+1)*n));
if ep == 0 & j == 1
    M1(1:n, 1:n) = coeff(L); % First Block M1
    M0(1:n, 1:n) = .5*coeff(L-1); % First Block M0
    M0(n+1:2*n, 1:n) = .5*(coeff(L-2)-2*coeff(L)); % Second Block M0
    for i = 3:k-1 % Blocks 3 to k-1 of M0
        M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = .5*(coeff(L-i)-coeff(L-i+2));
    end
    M0((k-1)*n+1:k*n, 1:n) = .5*(2*coeff(1)-coeff(3));
elseif ep == 0 & j == 2
    M1(1:n, 1:n) = 2*coeff(L); % First Block M1
    M0(1:n, 1:n) = coeff(L-1); % First Block M0
    M0(n+1:2*n, 1:n) = coeff(L-2)-coeff(L); % Second Block M0
    for i = 3:k % Blocks 3 to k of M0
        M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
    end
elseif ep == 1
    M1(1:n, 1:n) = 2*coeff(L); % Top Left Block M1
    M0(1:2*n, 1:n) = [coeff(L-1); coeff(L-2)-coeff(L)]; % First Two Left Block Rows M0
    for i = 3:ep-1 % Left Block Rows 3 to k - ep M0 
        M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
    end
    for i = 1:ep-1 % Right Block Rows 1 to k - ep - 1 M0
        M0((i-1)*n+1:i*n, 1:n) = -coeff(L-i-1);
    end
    
elseif ep >= 2 & ep <= k - 3
    
elseif ep == k - 2
    
else % ep == k - 1
    
end
    
