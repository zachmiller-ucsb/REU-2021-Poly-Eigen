function [data, Polys] = lshaped_varyepsilon(k, j, t) 
    %Tests effect of varying epsilons on K(L)/K(P) by generating t random
    %polynomials of degree k in j = type 1 or 2
    tic
    data = cell(1,t); %this is where the ratios of condition numbers for each eigenvalue are stored
    Polys = cell(1, t); %this is where the polynomials are stored
    
    for i = 1:t %Generate data for t polynomials
        
    P = create(k, j); %Makes a polynomial of degree k, type j.
    Polys{i} = P; %Stores P in Polys
    y = zeros(1,k); %Initializes y, a vector of size k
    
    %% for each epsilon between 0 and k-1
    for e = 0:(k-1)  
        [L1,L0] = lshaped(P,k,e,j); %Fill in maximum condition numbers
        collect = cnumber(L1,L0,k,e,j); %column vector filled with ratios of condition numbers for each eigenvalue for the given epsi
        y(1,e+1) = max(collect); %row vector with the maximum condition ratio for the given epsilon
    end  
    data{i} = y;
    fprintf("Finished a polynomial!\n")
    end
    toc
    
end

