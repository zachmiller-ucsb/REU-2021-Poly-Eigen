function [C1,C0] = C(n,a,j) %Constructing the C_epsilon and C_mu
%Construction of the minimal basis for block-minimal bases Chebyshev
%pencils.
% Inputs: n=size of blocs, j=1 or 2 (for Chebyshev of type 1 or 2),
% a=number of rows of the minimal basis.

        C1 = sym(zeros(a*n,(a+1)*n));
        C0 = sym(zeros(a*n,(a+1)*n));
        
        if a >= 1
        for i=0:a-1
            C0(i*n+1:(i+1)*n,i*n+1:(i+1)*n) = eye(n);
            C1(i*n+1:(i+1)*n,(i+1)*n+1:(i+2)*n) = -2*eye(n);
        end
        
        C1((a-1)*n+1:a*n,a*n+1:(a+1)*n) = -j*eye(n);
 
        for i=0:a-2
            C0(i*n+1:(i+1)*n,(i+2)*n+1:(i+3)*n) = eye(n);
        end
        end
end
