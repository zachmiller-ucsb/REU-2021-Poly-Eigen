function [L1,L0] = cPencil(M1,M0,j)
%Given M for a matrix polynomial and j = 1,2 specifying type 1 or 2 
%creates a Chebyshev pencil linearization
 
%specifies size of the coefficients (n = 1 means scalar coefficients)
n = 1; 

%Obtains mu and epsilon
[a,b] = size(M1); 
m = a/n-1;
e = b/n-1; 

%Defining L1 and L0
L1 = sym(zeros((m+e+1)*n)); 
L0 = sym(zeros((m+e+1)*n));
 
%Obtaining C_epsilon and C_mu
[Ce1,Ce0] = C(n, e,j);
[Cm1,Cm0] = C(n, m,2);
 
%Constructing the linearization
L1(1:(m+1)*n,1:(e+1)*n) = M1; 
L1((m+1)*n+1:(m+e+1)*n,1:(e+1)*n) = Ce1;
L1(1:(m+1)*n,(e+1)*n+1:(m+e+1)*n) = transpose(Cm1);
 
L0(1:(m+1)*n,1:(e+1)*n) = M0;
L0((m+1)*n+1:(m+e+1)*n,1:(e+1)*n) = Ce0;
L0(1:(m+1)*n,(e+1)*n+1:(m+e+1)*n) = transpose(Cm0);
 
end



