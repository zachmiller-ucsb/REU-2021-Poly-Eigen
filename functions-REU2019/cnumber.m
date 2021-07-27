function [cnumbers] = cnumber(L1, L0, k, e, j)
%Computes condition number ratio of linearization 
%to polynomial given epsilon, j = type 1,2
digits 40;
 
L0 = vpa(L0);
L1 = vpa(L1);
 
[V,~] = eig(-L1\L0); %columns of V are right eigenvectors, 
%diagonal entries of D are eigenvalues
 
[W,D] = eig((-L0/L1)'); %columns of W are left eigenvectors %%% we should use complex conjugate here in case we work with complex numbers
 
%k = length(P) - 1; %Obtains degree of polynomial
q = length(D); % This is the number of eigenvalues/eigenvectors (since they are singular) 
m = k-e-1; %Obtains mu
 
Vp = sym(zeros(1,q)); %Initialize recovery of eigenvectors
Wp = sym(zeros(1,q));
 
for i = 1:q %Recovery of right eigenvectors
    Vp(i) = V(e+1,i);  
end
 %i = 1:length(Vp)
 
for i = 1:q %Recovery of left eigenvectors
    Wp(i) = W(m+1,i);
end
 
cnumbers = sym(zeros(q,1)); %Set up vector
%that will contain condition numbers
 
for i=1:q %Fill in condition numbers
   ev = D(i,i); %eigenvalue
   
   Hx = V(:,i); %right and left eigenvectors of L
   yG = W(:,i);
   
   x = Vp(:,i); %right and left eigenvectors of P
   y = Wp(:,i);
    
   cheby = zeros(1,k+1); %Sum of Chebyshev polynomials 
   maxnorm = max(norm(L1),norm(L0));
   if j == 1
       for p=0:k
       cheby(1,p+1) = norm(chebyshevT(p,ev));
       end
   end
   
   if j == 2
       for p=0:k
       cheby(1,p+1) = norm(chebyshevU(p,ev));
       end  
   end
   
   l = j*abs(ev) + 1;   
   
   kl = maxnorm*l*norm(Hx)*norm(yG);
   kp = sum(cheby)*norm(x)*norm(y);
   
   cnumbers(i,1) = vpa(kl/kp);
end
end 



