% File that will run experiments for chebyshev backward error

d = 4; % degree
k = d;
polysize = 2; % size of polynomial
n = polysize;
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

j = 1; % Chebyshev Type 1 or 2
ep = 0; % Value of Epsilon for Linearization

%% Generate Polynomial in Monomial Basis

[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);



%% Perform Basis Conversion

Pcheb = sym(zeros(size(Pmon))); % This works without zero coefficients FIXED!
for r = 1:polysize
    for c = 1:polysize
        monEntry = coeffs(Pmon(r, c),'All'); % Returns coefficients with leading coefficient first 
        
        chebEntryPoly = mon2cheb(monEntry, j); % Returns coefficients with leading coefficient first in row vector 
        
        chebEntry = poly2sym(chebEntryPoly); % Returns symbolic polynomial taking leading coefficient first in row vector
        Pcheb(r, c) = chebEntry;
    end
end

coeff = polygen_split_smith_mtx(d, Pcheb); % ith entry is P_{i-1}


%% SCALING THE POLYNOMIAL
%Here, we scale P so that max\{norm(P4),...,norm(P0)\}=1

% For what follows
k = d;
n = polysize;

coeff = double(coeff); % Calculations no longer symbolic
coeffscal=zeros(n,n,d+1);
norms=zeros(k+1,1);
for i=1:k+1
    norms(i,1)=norm(coeff(:,:,i));
end
nmax=max(norms);
for i=1:k+1
    coeffscal(:,:,i)=coeff(:,:,i)/nmax;
end


%% Constructing the Linearization

[M1,M0] = Msubfamily(d,polysize,coeffscal,ep,j);
    
M1 = block2notblock(M1); 
M0 = block2notblock(M0);


[C1,C0] = cPencil(M1,M0,j,polysize,ep,d);

%% EIGENVALUE/EIGENVECTOR COMPUTATIONS
disp('Computing eigenvalues')

[Vc,ec]=eig(C0,-C1);
[ec,ind] = sort(diag(ec),'ascend');
Vc = Vc(:,ind);

% norm(sort(ec,'ascend') - sort(evs,'ascend'))/norm(evs)

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
    numC = norm((C0+C1*ec(i))*Vc(:,i));
    denC = norm(Vc(:,i))*max([nC0 nC1])*(j*abs(ec(i))+1);
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
        chebs = chebyshevT([0:d], ec(i)); % Generates vector of chebyshev polynomials 0:d evaluated at ec(i) (ith eigenvalue)
    else
        chebs = chebyshevU([0:d], ec(i));
    end
    resc=zeros(n,n);
    for j=1:d+1
        resc=resc+coeffscal(:,:,j)*chebs(j);
    end
    rc = norm(resc*Xc(:,i)); % Remainder Term

    atilde = sum(abs(chebs)); % Alpha Tilde Term
    
    back_error_Pc(i) = rc/(atilde*norm(Xc(:,i))); % ith entry in column vector of back error
    
    vector_norm_ratioc(i) = norm(Vc(:,i))/norm(Xc(:,i)); % ith entry in column vector of norm ratios
end

back_error_Pc

%% PLOTS

% BACKWARD ERROR RATIOS

figure
semilogy(back_error_Pc./back_error_C,'rx')

hold on

semilogy(vector_norm_ratioc,'bo')

hold on

legend('Pc/C')

title('Degree d - backward')

