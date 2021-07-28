function B = block2notblock(A)

% Function converts block matrix A to matrix B which is not in block form
% (for purposes of matlab computations)

sizeTemp = size(A);
totBlockRows = sizeTemp(3); % Total number block rows in A
totBlockCols = sizeTemp(4); % '' '' '' columns '' ''
m = sizeTemp(1); % Number rows in a block entry in A
n = sizeTemp(2); % '' columns '' '' '' '' '' ''
rows = sizeTemp(1)*sizeTemp(3); % Number rows in A
cols = sizeTemp(2)*sizeTemp(4); % Number columns in A
B = zeros(rows,cols); % Create B

for R = 1:totBlockRows
    for C = 1:totBlockCols
        B((R-1)*m+1:R*m, (C-1)*n+1:C*n) = A(:,:,R,C);
    end
end