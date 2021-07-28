% File that will run experiments for chebyshev backward error

d = 4; % degree
polysize = 3; % size of polynomial
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

symbolic = 1;

j = 1; % Chebyshev Type 1 or 2

ep = 2;

[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);

% test = polygen_split_smith_mtx(d, Pmon);
% test = double(test);
% P0 = test(:,:,1);
% P1 = test(:,:,2);
% P2 = test(:,:,3);
% P3 = test(:,:,4);
% P4 = test(:,:,5);
% polyeig(P0,P1,P2,P3,P4)
% evs

Pcheb = sym(zeros(size(Pmon)));
for r = 1:polysize
    for c = 1:polysize
        monEntry = coeffs(Pmon(r, c)); % Returns coefficients with leading coefficient last
        monEntry = fliplr(monEntry); % Returns coefficients with leading coefficient first (as mon2cheb takes it)
        chebEntryPoly = transpose(mon2cheb(monEntry, j)); % Returns coefficients with leading coefficient first IN ROW VECTOR
        chebEntryPoly = fliplr(chebEntryPoly); % Returns coefficients with leading coefficient last
        chebEntry = poly2sym(chebEntryPoly);
        Pcheb(r, c) = chebEntry;
    end
end

coeff = polygen_split_smith_mtx(d, Pcheb);

if symbolic == 1
    [M1sym,M0sym] = Msubfamily(d,polysize,coeff,ep,j);

    M1sym = block2notblock(M1sym);
    M0sym = block2notblock(M0sym);

    [C1sym,C0sym] = cPencil(M1sym,M0sym,j,polysize,ep,d);
    polyeig(C0sym,C1sym)
    evs
else
    coeff = double(coeff);

    [M1,M0] = Msubfamily(d,polysize,coeff,ep,j);
    M1 = block2notblock(M1);
    M0 = block2notblock(M0);

    [C1,C0] = cPencil(M1,M0,j,polysize,ep,d);
    
end

%% SCALING THE POLYNOMIAL
%Here, we scale P so that max\{norm(P4),...,norm(P0)\}=1

% For what follows
k = d;
n = polysize;

coeffscal=zeros(n,n,d+1);
nvect=zeros(n,1);
for i=1:(k+1)
    nvect(i,1)=norm(coeff(:,:,i));
end
nmax=max(nvect);
for i=1:k+1
    coeffscal(:,:,i)=coeff(:,:,i)/nmax;
end

%% EIGENVALUE/EIGENVECTOR COMPUTATIONS
disp('Computing eigenvalues')

[Vc,ec]=eig(C0,C1);
[ec,ind] = sort(diag(ec),'ascend');
Vc = Vc(:,ind); 

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

Xc = zeros(n,d*n);

for i=1:d*n
    Xc(:,i) = Vc(ep*n+1:(ep+1)*n,i);
end

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
        resc=resc+coeff(:,:,j)*chebs(j);
    end
    rc = norm(resc*Xc(:,i)); % Remainder Term
    
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

legend('Pc/C')

title('Degree d - backward')

    
   


