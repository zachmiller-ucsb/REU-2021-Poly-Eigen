function [M1,M0] = Msubfamily(d, polysize,coeff,ep,j)

% Constructs body of subfamily linearization for a matrix polynomial P
% written in the chebyshev basis dependent on parameters 0 <= ep <= k - 1 (size), k
% (degree of P), and j for chebyshev basis (type 1 --> j = 1, type 2 --> j = 2)

% Matrix coefficients A_i of P stored in block vector coeff(:,:, i + 1) bc
% 1 indexed

n = polysize; % Size of P 
sizeTemp = size(coeff); % Block vector size = [rows per block, columns per block, number of blocks]
L = sizeTemp(3); % NUMBER OF BLOCKS
k = d;
%M1 = zeros((k-ep)*n, (ep+1)*n);
M1 = zeros(n, n, k-ep, ep+1);
%M0 = zeros((k-ep)*n, (ep+1)*n);
M0 = zeros(n, n, k-ep, ep+1);
if ep == 0 & j == 1
    % M1(1:n, 1:n) = coeff(:,:,L); % First Block M1
    M1(:,:, 1) = coeff(:, :, L);
    %%%%
    % M0(1:n, 1:n) = .5*coeff(L-1); % First Block M0
    M0(:,:, 1) = .5*coeff(:,:, L-1);
    %%%
    %M0(n+1:2*n, 1:n) = .5*(coeff(L-2)-2*coeff(L)); % Second Block M0
    M0(:,:, 2) = .5*(coeff(:,:, L-2)-2*coeff(:,:, L));
    %%%
    for i = 3:k-1 % Blocks 3 to k-1 of M0
        %M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = .5*(coeff(L-i)-coeff(L-i+2));
        M0(:,:, i) = .5*(coeff(:,:, L-i)-coeff(:,:, L-i+2));
    end
    %M0((k-1)*n+1:k*n, 1:n) = .5*(2*coeff(1)-coeff(3));
    M0(:,:, k) = .5*(2*coeff(:,:, 1)-coeff(:,:, 3));
    
elseif ep == 0 & j == 2
    % M1(1:n, 1:n) = 2*coeff(L); % First Block M1
    M1(:, :, 1) = 2*coeff(:,:, L);
    
    % M0(1:n, 1:n) = coeff(L-1); % First Block M0
    M0(:, :, 1) = coeff(:,:, L-1);
    
    % M0(n+1:2*n, 1:n) = coeff(L-2)-coeff(L); % Second Block M0
    M0(:, :, 2) = coeff(:,:, L-2)-coeff(L);
    
    for i = 3:k % Blocks 3 to k of M0
        % M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
        M0(:, :, i) = coeff(:,:, L-i);
    end
    
elseif ep == 1
    % M1(1:n, 1:n) = 2*coeff(L); % Top Left Block M1
    M1(:, :, 1, 1) = 2*coeff(:,:, L);
    
    % M0(1:2*n, 1:n) = [coeff(L-1); coeff(L-2)-coeff(L)]; % First Two Left Block Rows M0
    M0(:, :, 1, 1) = coeff(:,:, L-1);
    M0(:, :, 2, 1) = coeff(:,:, L-2)-coeff(:,:, L)
    
    for i = 3:ep-1 % Left Block Rows 3 to k - ep M0 
        % M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
        M0(:, :, i, 1) = coeff(:,:, L-i);
    end
    for i = 1:k-ep-1 % Right Block Rows 1 to k-ep-1 M0
        % M0((i-1)*n+1:i*n, n+1:2*n) = -coeff(L-i+1);
        M0(:, :, i, 2) = -coeff(:,:, L-i+1);
    end
    % M0((k-ep-1)*n+1:(k-ep)*n, 2*n) = coeff(ep) - coeff(ep+2); % Right Block Row k - ep
    M0(:, :, k-ep, 2) = coeff(:,:, ep) - coeff(:,:, ep+2);
    
elseif ep >= 2 & ep <= k - 3
    % First Two columns same as ep == 1 in M0 and M1
    % M1(1:n, 1:n) = 2*coeff(L); 
    M1(:, :, 1, 1) = 2*coeff(:,:,L);
    
    % M0(1:2*n, 1:n) = [coeff(L-1); coeff(L-2)-coeff(L)]; 
    M0(:, :, 1, 1) = coeff(:,:,L-1);
    M0(:, :, 2, 1) = coeff(:,:,L-2)-coeff(:,:,L);
    
    for i = 3:ep-1 % Left Block Rows 3 to k - ep M0 
        % M0(2*n+(i-3)*n+1:2*n+(i-2)*n, 1:n) = coeff(L-i);
        M0(:, :, i, 1) = coeff(:,:, L-i);
    end
    for i = 1:k-ep-1 % Right Block Rows 1 to k-ep-1 M0
        % M0((i-1)*n+1:i*n, n+1:2*n) = -coeff(L-i+1);
        M0(:, :, i, 2) = -coeff(:,:, L-i+1);
    end
    
    % M0((k-ep-1)*n+1:(k-ep)*n, 2*n) = coeff(ep) - coeff(ep+2); 
    M0(:, :, k-ep, 2) = coeff(:,:, ep) - coeff(:,:, ep+2);
    % First two columns the same as ep == 1 in M0 and M1
    
    for i = 3:ep+1 % Last row, columns 3 to ep + 1
        % M0((k-ep-1)*n+1:(k-ep)*n, (i-1)*n+1:i*n) = coeff(ep+2-i);
        M0(:,:, k-ep, i) = coeff(:,:,ep+2-i);
    end
    
elseif ep == k - 2
    % M1(1:n, 1:n) = 2*coeff(L); % Top left block M1
    M1(:, :, 1, 1) = 2*coeff(:, :, L);
    
    % First Row M0
    % M0(1:n, 1:2*n) = [coeff(L-1) -coeff(L)]; % First Two Entries (Rest are Zeros)
    M0(:, :, 1, 1) = coeff(:,:,L-1);
    M0(:, :, 1, 2) = -coeff(:,:,L);
    
    % Second (Last) Row M0
    % M0(n+1:2*n, 1:2*n) = [coeff(ep+1) coeff(ep)-coeff(ep+2)]; % First Two Entries
    M0(:, :, 2, 1) = coeff(:,:,ep+1);
    M0(:, :, 2, 2) = coeff(:,:,ep)-coeff(:,:,ep+2);
    
    for i = 3:ep+1 % Columns 3 to ep + 1 (same as previous case)
        % M0((k-ep-1)*n+1:(k-ep)*n, (i-1)*n+1:i*n) = coeff(ep+2-i);
        M0(:, :, k-ep, i) = coeff(:,:,ep+2-i);
    end
    
else % ep == k - 1
    % M1(1:n, 1:n) = 2*coeff(L); % First Entry M1
    M1(:, :, 1, 1) = 2*coeff(:,:,L);
    
    % M0(1:n, 1:2*n) = [coeff(L-1) coeff(L-2)-coeff(L)]; % First Two Entries M0
    M0(:, :, 1) = coeff(:,:,L-1);
    M0(:, :, 2) = coeff(:,:,L-2)-coeff(:,:,L);
    
    for i = 3:k % Entries 3 to k in M0
        % M0(1:n, (i-1)*n+1:i*n) = coeff(L-i);
        M0(:, :, i) = coeff(:,:,L-i);
    end
end
    
