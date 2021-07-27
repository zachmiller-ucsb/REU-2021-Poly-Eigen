function [L1,L0] = lshaped(P, k, e, j)
%Creates L-shaped linearization for given P, epsilon, and mu. Notice that
%m=0 makes this the colleague pencil. j = type 1,2
 
 
m = k-e-1; %mu
 
M0 = sym(zeros(m+1,e+1)); %Set up M1 and M0
M1 = sym(zeros(m+1,e+1));
 
M1(1,1) = 2*P(1,1); %Start of algorithm
 
if m >= 1
    M0(2,1) = M0(2,1) - P(1,1);
end
 
if e >= 1
    M0(1,2) = M0(1,2) - P(1,1);
end
 
for i = 1:m
    M0(m-i+1,1) = M0(m-i+1,1) + P(1,k+1-e-i);
    if e >= 1
    M0(m-i+2,2) = M0(m-i+2,2) - P(1,k+1-e-i);
    end
end
 
for i = 0:e
    M0(m+1,e-i+1) = M0(m+1,e-i+1) + P(1,k+1-i);
end %End of algorithm
 
[L1,L0] = cPencil(M1,M0,j); %Puts M into a Chebyshev pencil linearization
 
end

