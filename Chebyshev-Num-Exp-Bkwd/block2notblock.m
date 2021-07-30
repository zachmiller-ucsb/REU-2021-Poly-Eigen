function B = block2notblock(A)

% Function converts block matrix A to matrix B which is not in block form
% (for purposes of matlab computations) 

% If block matrix has one dimension automatically combines to column vector

sizeTemp = size(A);
metaSize = size(size(A)); % Want to check if A is block row or column vector
m = sizeTemp(1); % Number rows in a block entry in A
n = sizeTemp(2); % '' columns '' '' '' '' '' ''
if metaSize(2) == 3
    totBlockCols = 1;
    totBlockRows = sizeTemp(3);
else
    totBlockRows = sizeTemp(3); % Total number block rows in A
    totBlockCols = sizeTemp(4); % '' '' '' columns '' ''
end
rows = sizeTemp(1)*totBlockRows; % Number rows in A
cols = sizeTemp(2)*totBlockCols; % Number columns in A
B = zeros(rows,cols); % Create B

for R = 1:totBlockRows
    for C = 1:totBlockCols
        B((R-1)*m+1:R*m, (C-1)*n+1:C*n) = A(:,:,R,C);
    end
end