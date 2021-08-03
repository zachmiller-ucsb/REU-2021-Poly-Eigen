% File that will run experiments for chebyshev backward error

d = 4; % degree
polysize = 2; % size of polynomial
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

symbolic = 0; % Doesn't work at the moment, because we double the matrix polynomial before constructing linearization

j = 1; % Chebyshev Type 1 or 2

ep = 3;

[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);

% test = polygen_split_smith_mtx(d, Pmon);
% test = double(test)
% P0 = test(:,:,1);
% P1 = test(:,:,2);
% P2 = test(:,:,3);
% P3 = test(:,:,4);
% P4 = test(:,:,5);
% Tevs = polyeig(P0,P1,P2,P3,P4)
% evs





Pcheb = sym(zeros(size(Pmon)));
for r = 1:polysize
    for c = 1:polysize
        monEntry = coeffs(Pmon(r, c)); % Returns coefficients with leading coefficient last % SUCCESS
        monEntry = fliplr(monEntry); % Returns coefficients with leading coefficient first (as mon2cheb takes it) % SUCCESS
        
        chebEntryPoly = mon2cheb(monEntry, j); % Returns coefficients with leading coefficient first in row vector % SUCCESS

%         chebEntryPoly = fliplr(chebEntryPoly); % Returns coefficients
%         with leading coefficient last % SUCCESS  --- This was the
%         error

        chebEntry = poly2sym(chebEntryPoly);
        Pcheb(r, c) = chebEntry;
    end
end

%%Another Test
% double(coeffs(Pmon(1,1)))
% double(coeffs(Pcheb(1,1)))

coeff = polygen_split_smith_mtx(d, Pcheb);


%%test using conversion back to Pmon from Pcheb
% Pmoncoeff = polygen_split_smith_mtx(d, Pmon);
% Pmoncoeff = double(Pmoncoeff)
% doubPmon = zeros(polysize);
% for i=1:d+1
%     doubPmon = doubPmon + Pmoncoeff(:,:,i)*sym('x');
% end
% doubPmon
% coeff = double(coeff);


% Check Change of Basis Degree 4 (Chebyshev Polynomials)
back2mon = zeros(polysize);
for i=1:d+1
    back2mon = back2mon + coeff(:,:,i)*chebyshevT(i-1,sym('x'));
end
back2mon = polygen_split_smith_mtx(d, back2mon);
back2mon = double(back2mon)
P0 = back2mon(:,:,1);
P1 = back2mon(:,:,2);
P2 = back2mon(:,:,3);
P3 = back2mon(:,:,4);
P4 = back2mon(:,:,5);
Back2Monevs = polyeig(P0,P1,P2,P3,P4);
norm(sort(evs,'ascend') - sort(Back2Monevs,'ascend'))/norm(evs)

%%

%% SCALING THE POLYNOMIAL
%Here, we scale P so that max\{norm(P4),...,norm(P0)\}=1

% For what follows
k = d;
n = polysize;

coeff = double(coeff);
coeffscal=zeros(n,n,d+1);
norms=zeros(k+1,1);
for i=1:(k+1)
    norms(i,1)=norm(coeff(:,:,i));
end
nmax=max(norms);
for i=1:k+1
    coeffscal(:,:,i)=coeff(:,:,i)/nmax;
end

%% BUILDING THE LINEARIZATION

if symbolic == 1 % This will not run
    [M1sym,M0sym] = Msubfamily(d,polysize,coeffscal,ep,j);

    M1sym = block2notblock(M1sym);
    M0sym = block2notblock(M0sym);

    [C1sym,C0sym] = cPencil(M1sym,M0sym,j,polysize,ep,d);
 
else
    coeffscal = double(coeffscal); % Check this

    [M1,M0] = Msubfamily(d,polysize,coeffscal,ep,j);
    
    if ep == d - 1
        M1 = pagetranspose(M1);
        M0 = pagetranspose(M0);
        M1 = block2notblock(M1);
        M0 = block2notblock(M0);
%         M1 = transpose(M1);
%         M0 = transpose(M0);
    else
        M1 = block2notblock(M1); % Fix this for ep = k - 1
        M0 = block2notblock(M0);
    end
    

    [C1,C0] = cPencil(M1,M0,j,polysize,ep,d);
    
end
%% EIGENVALUE/EIGENVECTOR COMPUTATIONS
disp('Computing eigenvalues')

[Vc,ec]=eig(C0,-C1);
[ec,ind] = sort(diag(ec),'ascend');
Vc = Vc(:,ind);
ec - sort(evs,'ascend')

% disp('Check Linearization:')
% ec
% evs
% sort(evs,'ascend')

%% Eigenvector Recovery

Xc = zeros(n,d*n);

for i=1:d*n
    Xc(:,i) = Vc(ep*n+1:(ep+1)*n,i);
end

%% LINEARIZED BACKWARD ERROR
disp('Linearizations backward errors')

nC0 = norm(C0);
nC1 = norm(C1);

%Compute backward errors
back_error_C = zeros(d*n,1);
for i=1:d*n
    numC = norm((C0*ec(i)+C1)*Vc(:,i));
    denC = norm(Vc(:,i))*max([nC0 nC1])*(abs(ec(i)+1));
    back_error_C(i) = numC/denC;
end

%% POLYNOMIAL BACKWARD ERRORS
disp('polynomial backward errors')

%BACKWARD ERRORS

back_error_Pc = zeros(d*n,1);

vector_norm_ratioc = zeros(d*n,1);

for i=1:d*n
    if r == 1
        chebs = chebyshevT([0:d], ec(i));
    else
        chebs = chebyshevU([0:d], ec(i));
    end
    resc=zeros(n,n);
    for j=1:d+1
        resc=resc+coeffscal(:,:,j)*chebs(j);
    end
    rc = norm(resc*Xc(:,i)); % Remainder Term
    
    % Check Eigenpairs Polynomial
    
%     disp('Check Eigenpairs Polynomial')
%     for i=1:d*polysize
%         resc*Xc(:,i)
%     end
    
    moduli = zeros(d+1,1);
    for j=1:d+1
        moduli(j,1) = abs(chebs(j));
    end
    atilde = sum(moduli); % Alpha Tilde Term
    
    back_error_Pc(i) = rc/(atilde*norm(Xc(:,i))); % ith entry in column vector of back error
    
    vector_norm_ratioc(i) = norm(Vc(:,i))/norm(Xc(:,i)); % ith entry in column vector of norm ratios
end
    
%% PLOTS

% BACKWARD ERROR RATIOS

figure
semilogy(back_error_Pc./back_error_C,'rx')

hold on

semilogy(vector_norm_ratioc,'bo')

hold on

% Plotting Bounds
% Backward Error Bound Polynomial
% if ep == 0 & j == 1
%     
% elseif ep == 0 & j == 2
%     
% elseif ep == k - 1 & j == 1
%     
% elseif ep == k - 1 & j == 2
%     
% else
%     
% end
% semilogy
% Ratio of EigenVector Bound
% semilogy()

legend('Pc/C')

title('Degree d - backward')
   


