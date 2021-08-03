% Test File for Linearization Construction 2 x 2 case, degree 5

%% Generates simple polynomial to test structure
testPoly = zeros(2,2,6);
for i=1:6
    testPoly(:,:,i) = [i-1 i-1;i-1 i-1];
end

%% Runs Body Generation

% ep = 0 j = 1
disp('ep = 0, j = 1')

[M1,M0] = Msubfamily(5,2,testPoly,0,1);
M1 = block2notblock(M1)
M0 = block2notblock(M0)

% ep = 0 j = 2
disp('ep = 0 j = 2')

[M1,M0] = Msubfamily(5,2,testPoly,0,2);
M1 = block2notblock(M1)
M0 = block2notblock(M0)

% ep = 1 
disp('ep = 1')

[M1,M0] = Msubfamily(5,2,testPoly,1,1);
M1 = block2notblock(M1)
M0 = block2notblock(M0)

% ep = 2
disp('ep = 2 (2<=ep<=k-3)')

[M1,M0] = Msubfamily(5,2,testPoly,2,1);
M1 = block2notblock(M1)
M0 = block2notblock(M0)

% ep = 3
disp('ep = 3 (k-2)')

[M1,M0] = Msubfamily(5,2,testPoly,3,1);
M1 = block2notblock(M1)
M0 = block2notblock(M0)

% ep = 4 
disp('ep = 4 (k-1)')

[M1,M0] = Msubfamily(5,2,testPoly,4,1);

M1 = block2notblock(M1)
M0 = block2notblock(M0)

%% Whole Pencil

[L1,L0] = cPencil(M1,M0,1,2,4,5);

L1
L0


% %% Test Change of Basis WORKS!!! (As far as I can tell)
% 
% % Generates simple polynomial to test structure
% testPoly = sym(zeros(2,2,10));
% for i=1:10
%     testPoly(:,:,i) = [i-1 i-1;i-1 i-1];
% end
% 
% Pmon = sym(zeros(2,2));
% for i=1:10
%     Pmon = Pmon + testPoly(:,:,i)*sym('x')^(i-1);
% end
% 
% Pcheb = sym(zeros(size(Pmon))); % This works without zero coefficients FIXED
% for r = 1:2
%     for c = 1:2
%         monEntry = coeffs(Pmon(r, c),'All'); % Returns coefficients with leading coefficient last 
%         
%         chebEntryPoly = mon2cheb(monEntry, 1); % Returns coefficients with leading coefficient first in row vector 
%         
%         chebEntry = poly2sym(chebEntryPoly); % Returns symbolic polynomial taking leading coefficient first in row vector
%         Pcheb(r, c) = chebEntry;
%     end
% end
% coeff = polygen_split_smith_mtx(9, Pcheb);
% 
% back2mon = sym(zeros(2));
% for i=1:10
%     back2mon = back2mon + coeff(:,:,i)*chebyshevT(i-1,sym('x'));
% end
% back2mon = polygen_split_smith_mtx(2, back2mon)


% %% Test Generating Function
% 
% d = 1; % degree
% k = d;
% polysize = 2; % size of polynomial
% n = polysize;
% pkmean = 0; % mean for pseudosmith
% pkwidth = 50; % std for pseudosmith
% 
% j = 1; % Chebyshev Type 1 or 2
% ep = 1; % Value of Epsilon for Linearization
% 
% [evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);
% 
% test = polygen_split_smith_mtx(d, Pmon);
% P0 = test(:,:,1);
% P1 = test(:,:,2);
% % P2 = test(:,:,3);
% % P3 = test(:,:,4);
% % P4 = test(:,:,5);
% Tevs = eig(P0,-P1);
% (norm(sort(evs,'ascend') - sort(Tevs,'ascend')))/norm(evs)
% 
% %% Test Backward Error Computation
% 
% r = 1;
% 
% % Generates simple polynomial to test structure
% testPoly = sym(zeros(2,2,10));
% for i=1:5
%     testPoly(:,:,i) = [i-1 i-1;i-1 i-1];
% end
% 
% if r == 1
%     chebs = chebyshevT([0:4], 2);
% else
%     chebs = chebyshevU([0:4], 2);
% end
% resc=zeros(2,2);
% for j=1:5
%     resc=resc+testPoly(:,:,j)*chebs(j)
% end
