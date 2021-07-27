function [P]=create(k,j)
 
%Given degree k, this program generates a random vector with real entries in [-1,1] and then a polynomial in the monomial basis whose
%roots are those random numbers. Then, it expresses that polynomial in the
%Chebyshev of type 1 or 2 basis depending on j = 1,2
 la=2*rand([k,1])-1;
 B=sym(poly(la));
 P=mon2cheb(B,j);
 a=max(abs(P));
 P=P/a;
 P=transpose(P)  
end
