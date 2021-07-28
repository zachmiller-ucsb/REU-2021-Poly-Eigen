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
    for i = 1:k-ep-1 % Right Block Rows 1 to k-ep-1 M0
        M0((i-1)*n+1:i*n, n+1:2*n) = -coeff(L-i+1);
    end
    M0((k-ep-1)*n+1:(k-ep)*n, 2*n) = coeff(ep) - coeff(ep+2); % Right Block Row k - ep
    
elseif ep >= 2 & ep <= k - 3
    % First Two columns same as ep == 1 in M0 and M1
    M1(1:n, 1:n) = 2*coeff(L); 
    M0(1:2*n, 1:n) = [coeff(L-1); coeff(L-2)-coeff(L)]; 
    for i = 3:ep-1  
        M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
    end
    for i = 1:k-ep-1 
        M0((i-1)*n+1:i*n, n+1:2*n) = -coeff(L-i+1);
    end
    M0((k-ep-1)*n+1:(k-ep)*n, 2*n) = coeff(ep) - coeff(ep+2); 
    % First two columns the same as ep == 1 in M0 and M1
    
    for i = 3:ep+1 % Last row, columns 3 to ep + 1
        M0((k-ep-1)*n+1:(k-ep)*n, (i-1)*n+1:i*n) = coeff(ep+2-i);
    
elseif ep == k - 2
    M1(1:n, 1:n) = 2*coeff(L); % Top left block M1
    % First Row M0
    M0(1:n, 1:2*n) = [coeff(L-1) -coeff(L)]; % First Two Entries (Rest are Zeros)
    % Second (Last) Row M0
    M0(n+1:2*n, 1:2*n) = [coeff(ep+1) coeff(ep)-coeff(ep+2)]; % First Two Entries
    for i = 3:ep+1 % Columns 3 to ep + 1 (same as previous case)
        M0((k-ep-1)*n+1:(k-ep)*n, (i-1)*n+1:i*n) = coeff(ep+2-i);
    
else % ep == k - 1
    M1(1:n, 1:n) = 2*coeff(L); % First Entry M1
    M0(1:n, 1:2*n) = [coeff(L-1) coeff(L-2)-coeff(L)]; % First Two Entries M0
    for i = 3:k % Entries 3 to k in M0
        M0(1:n, (i-1)*n+1:i*n) = coeff(L-i);
end
    