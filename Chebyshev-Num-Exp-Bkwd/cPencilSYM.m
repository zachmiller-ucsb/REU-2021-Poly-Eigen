function [L1,L0] = cPencil(M1,M0,j,polysize,ep,d)
%Given M for a matrix polynomial and j = 1,2 specifying type 1 or 2 
%creates a Chebyshev pencil linearization

%specifies size of the coefficients (n = 1 means scalar coefficients)
n = polysize; 

% Obtains mu 
mu = d - 1 - ep;

%Defining L1 and L0
L1 = sym(zeros(d*n)); 
L0 = sym(zeros(d*n));
 
%Obtaining C_epsilon and C_mu
[Ce1,Ce0] = C(n, ep,j);
[Cm1,Cm0] = C(n, mu,2);
 
%Constructing the linearization
L1(1:(mu+1)*n,1:(ep+1)*n) = M1; 
L1((mu+1)*n+1:d*n,1:(ep+1)*n) = Ce1;
L1(1:(mu+1)*n,(ep+1)*n+1:d*n) = transpose(Cm1);
 
L0(1:(mu+1)*n,1:(ep+1)*n) = M0;
L0((mu+1)*n+1:d*n,1:(ep+1)*n) = Ce0;
L0(1:(mu+1)*n,(ep+1)*n+1:d*n) = transpose(Cm0);
 
end



