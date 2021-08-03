% File that will run experiments for chebyshev backward error

d = 15; % degree
k = d;
polysize = 2; % size of polynomial
n = polysize;
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

j = 1; % Chebyshev Type 1 or 2
ep = 0; % Value of Epsilon for Linearization

%% Constructing Sample Polynomial 
% Generate Matrix Polynomial Pmon in Monomial Basis with Eigenvalues evs in [-1,1]
[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth); 

% Generate Change of Basis Matrix from Chebyshev Basis to Monomial Basis
L = std_xchg_chebyshev_monomial(d,j); 

% Obtain Matrix Coefficients of Pmon
PmonCoeff = polygen_split_smith_mtx(d, Pmon);

% Perform Change of Basis on Pmon's Matrix Coefficients
Pcheb = polygen_rewrite_polynomial(L, PmonCoeff);

Pcheb2Mon = sym(zeros(n));
for i=1:d+1
    Pcheb2Mon = Pcheb2Mon + Pcheb(:,:,i)*chebyshevT(i-1,sym('x'));
end


%% Scaling the Polynomial
%Here, we scale P so that max\{norm(Pd),...,norm(P0)\}=1

Pcheb = double(Pcheb);
PchebScal=zeros(n,n,d+1);
norms=zeros(d+1,1);
for i=1:(d+1)
    norms(i,1)=norm(Pcheb(:,:,i));
end
nmax=max(norms);
for i=1:d+1
    PchebScal(:,:,i)=Pcheb(:,:,i)/nmax;
end

%% Constructing the Linearization

[M1,M0] = Msubfamily(d,polysize,PchebScal,ep,j);
    
if ep == d - 1 % Fix this for ep = d - 1
    M1 = pagetranspose(M1);
    M0 = pagetranspose(M0);
    M1 = block2notblock(M1);
    M0 = block2notblock(M0);
%         M1 = transpose(M1);
%         M0 = transpose(M0);
else
    M1 = block2notblock(M1); 
    M0 = block2notblock(M0);
end

[C1,C0] = cPencil(M1,M0,j,polysize,ep,d);

%% EIGENVALUE/EIGENVECTOR COMPUTATIONS
disp('Computing eigenvalues')

[Vc,ec]=eig(C0,-C1);
[ec,ind] = sort(diag(ec),'ascend');
Vc = Vc(:,ind);

ec - sort(evs,'ascend')
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
    r = j;
    if r == 1
        chebs = chebyshevT([0:d], ec(i));
    else
        chebs = chebyshevU([0:d], ec(i));
    end
    resc=zeros(n,n);
    for j=1:d+1
        resc=resc+PchebScal(:,:,j)*chebs(j);
    end
    rc = norm(resc*Xc(:,i)); % Remainder Term
    % resc*Xc(:,i)
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
