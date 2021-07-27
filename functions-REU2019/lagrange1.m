function  lagrange(k, n, )

%Uses algorithm 1 to find the eigenvalues of nonlinear eigenvalue problem
%T(lambda)x=0 called "Gun" in the nlevp collection. This function
%interpolates the function T at nodes  in the interval [11880, 11882] given
%by option n and then constructs the colleague pencil associated with the
%interpolating matrix polynomial P of degree k using Lagrange polynomials. 
% Finally, it computes the eigenvalues and eigenvectors of P in that
% interval and finds its backward error. 


[coeffs, FUN]=nlevp('gun'); % coeffs gives the matrix coefficients of the split form of T, that is, 
% T(la)= f1(la)C1 + f2(la)C2+ ... +fn(la)Cn. They are 9956x 9956 matrices.



x=zeros(k+1, 1);  % vector of interpolating nodes.
a1=11880+11882;
a2=2;
if n==0
    x(i)= a1/2+a2/2*cos((2i-1)/2*(k+1)*pi);  %scaled Tchebyshev nodes.
else
    x=11881+2*randn(k+1,1);                  % random nodes with normal distribution
end



Tx=cell(1,k+1); % evaluation of interpolating nodes at T. This is also the set of matrix coefficients for the 
                       %interpolating polynomial in Lagrange form.
for j=1:k+1
Tx{j}=(nlevp('eval', 'gun', x(j) )); 
end
%construction of the linearization

[L1, L0]= colleague(Tx, x); %construction of the linearization

[V, D]=eig(-L0,L1);  % D is a diagonal matrix containing the eigenvalues and the columns of V are the corresponding eigenvectors.
D=diag(D);
[~, I] = sort(abs(D));  % order eigenvalues in increasing order of modulus
    D = D(I);
    V1 = V1(:,I');      % order eigenvectors to match the order of the eigenvalues.


[r,t]=size(D);


berror=zeros(r,1);  % vector of backward errors. 
normcoeffs=zeros(4,1);
for j=1:4
    normcoeffs(j)=norm(coeffs{j});
end
alfa=max(normcoeffs);  %%% weight for error backward.

f=zeros(r, 1);
for j=1:r
f(j)=FUN(D(j));
f(j)=sum(f(j));
end

for j=1:r
    
   berror(i)=norm(nlevp('eval',gun,D(j))*V1(:,j))/(norm(V1(:,j))*alfa*f(j) ) ;
end
   
    
    


    