wstring = 'norm'; % 'norm' or 'coeff' %coeff:relative-relative; norm: relative-absolute
scaled = 'bad'; % 'well' or 'bad' %valor de rho
symbolic = 'no'; % 'no' or 'yes'
epsilon = 3; % 0 <= epsilon <= degree
r = 2; % Chebyshev... First Kind: 1 Second Kind: 2
%% DEGREE-4 RANDOM MATRIX POLYNOMIAL: Experiment for the paper
disp('Building matrix polynomial')
n = 20; %size of the polynomial

if strcmp(scaled,'bad')==1 
    %randn('state',1010)     % siempre mismo experimento
    P4 = 1e1*(randn(n)+sqrt(-1)*randn(n)); %1e1
    P3 = 1e6*(randn(n)+sqrt(-1)*randn(n));
    P2 = 1e7*(randn(n)+sqrt(-1)*randn(n));
    P1 = 1e6*(randn(n)+sqrt(-1)*randn(n));
    P0 = 1e4*(randn(n)+sqrt(-1)*randn(n)); %1e4
    

else  %%% rho =1
    %randn('state',1010)
    P4 = randn(n)+sqrt(-1)*randn(n);
    P3 = randn(n)+sqrt(-1)*randn(n);
    P2 = randn(n)+sqrt(-1)*randn(n);
    P1 = randn(n)+sqrt(-1)*randn(n);
    P0 = randn(n)+sqrt(-1)*randn(n);    
    
    [U,S,V]=svd(P4);
    S(end-5,end-5) = 1e-5;
    S(end-4,end-4) = 1e-5;
    S(end-3,end-3) = 1e-5;
    S(end-2,end-2) = 1e-5;
    S(end-1,end-1) = 1e-6;
    S(end,end)=1e-6;
    P4 = U*S*V';
    
    [U,S,V]=svd(P0);
    S(end-5,end-5) = 1e-5;
    S(end-4,end-4) = 1e-5;
    S(end-3,end-3) = 1e-5;
    S(end-2,end-2) = 1e-5;
    S(end-1,end-1) = 1e-6;
    S(end,end)=1e-6;
    P0 = U*S*V';


end



%% SCALING THE POLYNOMIAL
%Here, we scale P so that max\{norm(P4),...,norm(P0)\}=1
disp('Scaling the polynomial')
n4 = norm(P4);
n3 = norm(P3);
n2 = norm(P2);
n1 = norm(P1);
n0 = norm(P0);
nmax = max([n0 n1 n2 n3 n4]);
P4 = P4/nmax;
P3 = P3/nmax;
P2 = P2/nmax;
P1 = P1/nmax;
P0 = P0/nmax;
n4 = norm(P4);
n3 = norm(P3);
n2 = norm(P2);
n1 = norm(P1);
n0 = norm(P0);
nmax = 1;
%% LINEARIZATION
disp('Building the linearizations')
o = zeros(n);
id = eye(n);
%Dk linearization
Ak = [o o -P4 o; o -P4 -P3 o; -P4 -P3 -P2 o;o o o P0];
Bk = [o o o P4;o o P4 P3; o P4 P3 P2; P4 P3 P2 P1];
nAk = norm(Ak);
nBk = norm(Bk);
%D1 linearization
A1 = [P3 P2 P1 P0; P2 P1 P0 o; P1 P0 o o; P0 o o o];
B1 = [P4 o o o; o -P2 -P1 -P0; o -P1 -P0 o; o -P0 o o];
nA1 = norm(A1);
nB1 = norm(B1);
%Companion linearization
A  = [P3 P2 P1 P0;-id o o o; o -id o o; o o -id o];
B = [P4 o o o; o id o o; o o id o; o o o id];
%Lawrence-Perez Subfamily Linearization
if epsilon == 0
    K22T = [id o o; o id o; id o id; o id o];
    if r == 1
        M1 = [P4; o; o; o];
        M2 = [.5*P3; .5*P2-P4; .5*(P1-P3); P0-.5*P2];
        K21T = [o o o; -2*id o o; o -2*id o; o o -id];
    else
        M1 = [2*P4; o; o; o];
        M2 = [P3; P2-P4; P1; P0];
        K21T = [o o o; -2*id o o; o -2*id o; o o -2*id];
    end
    C1 = [M1 K21T];
    C2 = [M2 K22T];
elseif epsilon == 1
    M1 = [2*P4 o; o o; o o];
    M2 = [P3 -P4; P2-P4 -P3; P1 -P0-P2];
    K12 = [id o]; 
    K22T = [id o; o id; id o]; 
    if r == 1
        K11 = [o -id];
        K21T = [o o; -2*id o; o -id];
    else
        K11 = [o -2*id];
        K21T = [o o; -2*id o; o -2*id];
    end
    C1 = [M1 K21T; K11 o o];
    C2 = [M2 K22T; K12 o o];
elseif epsilon == 2
    M1 = [2*P4 o o; o o o];
    M2 = [P3 -P4 o; P2 P1-P3 P0];
    K12 = [id o id; o id o];
    K22T = [id; o];
    if r == 1
        K11 = [o -2*id o; o o -id];
        K21T = [o; -id];
    else
        K11 = [o -2*id o; o o -2*id];
        K21T = [o; -2*id];
    end
    C1 = [M1 K21T; K11 zeros(2*n,n)];
    C2 = [M2 K22T; K12 zeros(2*n,n)];
else
    M1 = [2*P4 o o o];
    M2 = [P3 P2-P4 P1 P0];
    K12 = [id o id o; o id o id; o o id o];
    if r == 1
        K11 = [o -2*id o o; o o -2*id o; o o o -id];
    else
        K11 = [o -2*id o o; o o -2*id o; o o o -2*id];
    end
    C1 = [M1; K11];
    C2 = [M2; K12];
end
nC1 = norm(C1);
nC2 = norm(C2);
%nA = norm(A);
%nB = norm(B);
%H linearization
%AH = [-P4 o o o; o P2 -i o; o -i o o; o o o P0];
%BH = [o P4 o o; P4 P3 o o; o o o i; o o i P1];
%nAH = norm(AH);
%nBH = norm(BH);
%G linearization
%AG = [o P0 o o; P0 P1 o o; o o o i; o o i P3];
%BG = [-P0 o o o; o P2 -i o; o -i o o; o o o P4];
%nAG = norm(AG);
%nBG = norm(BG);
%% EIGENVALUE/EIGENVECTOR COMPUTATIONS
disp('Computing eigenvalues')

[Vc,ec]=eig(C1,C2);
[ec,ind] = sort(diag(ec),'ascend');
ec
Vc = Vc(:,ind); 


% [Vk,ek]=eig(Ak,-Bk);
% [ek,ind] = sort(diag(ek),'ascend');
% Vk = Vk(:,ind);


% [V1,e1]=eig(A1,-B1);
% [e1,ind] = sort(diag(e1),'ascend');
% V1 = V1(:,ind);

%[VH,eH]=eig(AH,-BH);
%[eH,ind] = sort(diag(eH),'descend');
%VH = VH(:,ind);

%[VG,eG]=eig(AG,-BG);
%[eG,ind] = sort(diag(eG),'descend');
%VG = VG(:,ind);


%% LINEARIZED BACKWARD ERROR
disp('Linearizations backward errors')

%Compute backward errors
back_error_C = zeros(4*n,1);
%back_error_Dk = zeros(4*n,1);
%back_error_D1 = zeros(4*n,1);
%back_error_DH = zeros(4*n,1);
%back_error_DG = zeros(4*n,1);
for i=1:4*n
    numC = norm((C1*ec(i)+C2)*Vc(:,i));
    %numk = norm((Bk*ek(i)+Ak)*Vk(:,i));
    %num1 = norm((B1*e1(i)+A1)*V1(:,i));
    %numH = norm((BH*eH(i)+AH)*VH(:,i));
    %numG = norm((BG*eG(i)+AG)*VG(:,i));
    if strcmp(wstring,'coeff')==1 
      %denk = norm(Vk(:,i))*(nBk*abs(ek(i))+nAk);
      %den1 = norm(V1(:,i))*(nB1*abs(e1(i))+nA1);
      %denH = norm(VH(:,i))*(nBH*abs(eH(i))+nAH);
      %denG = norm(VG(:,i))*(nBG*abs(eG(i))+nAG);
    else
       denC = norm(Vc(:,i))*max([nC1 nC2])*(abs(ec(i)+1));
       % denk = norm(Vk(:,i))*max([nBk nAk])*(abs(ek(i))+1);
       % den1 = norm(V1(:,i))*max([nB1 nA1])*(abs(e1(i))+1);
       % denH = norm(VH(:,i))*max([nBH nAH])*(abs(eH(i))+1);
       % denG = norm(VG(:,i))*max([nBG nAG])*(abs(eG(i))+1);
    end
    back_error_C(i) = numC/denC;
    %back_error_Dk(i) = numk/denk;
	%back_error_D1(i) = num1/den1;
    %back_error_DH(i) = numH/denH;
    %back_error_DG(i) = numG/denG;
end

%% POLYNOMIAL BACKWARD ERRORS
disp('polynomial backward errors')
Xc = zeros(n,4*n);
%Xk = zeros(n,4*n);
%hk = zeros(1,4*n);
%X1 = zeros(n,4*n);
%%XH = zeros(n,4*n);
%XG = zeros(n,4*n);
%recovered eigenvectors
for i=1:4*n
    Xc(:,i) = Vc(epsilon*n+1:(epsilon+1)*n,i); %CHECK THIS 2021
    %Xk(:,i) = Vk(3*n+1:4*n,i); %recover from last block
    %X1(:,i) = V1(1:n,i); %recover from first block
   
   % if abs(eH(i))<=1
   %     XH(:,i) = VH(3*n+1:4*n,i); %recover from last block
    %else
   %     XH(:,i) = VH(n+1:2*n,i); %recover from second block
   % end
   % if abs(eG(i))>1
   %     XG(:,i) = VG(3*n+1:4*n,i); %recover from last block
   % else
    %    XG(:,i) = VG(n+1:2*n,i); %recover from second block
    %end
end

%BACKWARD ERRORS

back_error_Pc = zeros(4*n,1);
%back_error_Pk = zeros(4*n,1);
%back_error_P1 = zeros(4*n,1);
%back_error_PG = zeros(4*n,1);
%back_error_PH = zeros(4*n,1);

vector_norm_ratioc = zeros(4*n,1);
%vector_norm_ratio1 = zeros(4*n,1);
%vector_norm_ratiok = zeros(4*n,1);
for i=1:4*n
    % C
    if r == 1
        chebs = chebyshevT([0, 1, 2, 3, 4], ec(i));
    else
        chebs = chebyshevU([0, 1, 2, 3, 4], ec(i));
    end
    x0 = chebs(1);
    x1 = chebs(2); 
    x2 = chebs(3); 
    x3 = chebs(4); 
    x4 = chebs(5); 
    resc = P4*x4+P3*x3+P2*x2+P1*x1+P0*x0;
    rc = norm(resc*Xc(:,i));
    % Dk
%     resk = P4*ek(i)+P3; 
%     resk = ek(i)*resk+P2;
%     resk = ek(i)*resk+P1;
%     resk = ek(i)*resk+P0;
%     rk = norm(resk*Xk(:,i));
    % D1
%     res1 = P4*e1(i)+P3; 
%     res1 = e1(i)*res1+P2;
%     res1 = e1(i)*res1+P1;
%     res1 = e1(i)*res1+P0;
%     r1 = norm(res1*X1(:,i));
    % H
   % resH = P4*eH(i)+P3; 
   % resH = eH(i)*resH+P2;
   % resH = eH(i)*resH+P1;
   % resH = eH(i)*resH+P0;
   % rH = norm(resH*XH(:,i));
    % G
   % resG = P4*eG(i)+P3; 
   % resG = eG(i)*resG+P2;
   % resG = eG(i)*resG+P1;
   % resG = eG(i)*resG+P0;
   % rG = norm(resG*XG(:,i));
    
   y0 = abs(x0);
   y1 = abs(x1);
   y2 = abs(x2);
   y3 = abs(x3);
   y4 = abs(x4);
   %eek = abs(ek(i));
   %ee1 = abs(e1(i));
   % eeH = abs(eH(i));
   % eeG = abs(eG(i));
%     if strcmp(wstring,'coeff')==1
%         w0 = n0;
%         w1 = n1;
%         w2 = n2;
%         w3 = n3;
%         w4 = n4;
%     else
%         w0 = nmax;
%         w1 = nmax;
%         w2 = nmax;
%         w3 = nmax;
%         w4 = nmax;
%     end
    
    % C
    polyc = sum([y0, y1, y2, y3, y4]);
    
    % Dk
%     polyk = w4*eek+w3;
%     polyk = eek*polyk+w2;
%     polyk = eek*polyk+w1;
%     polyk = eek*polyk+w0;
    % D1
%     poly1 = w4*ee1+w3;
%     poly1 = ee1*poly1+w2;
%     poly1 = ee1*poly1+w1;
%     poly1 = ee1*poly1+w0;
    % G
   % polyG = w4*eeG+w3;
   % polyG = eeG*polyG+w2;
   % polyG = eeG*polyG+w1;
   % polyG = eeG*polyG+w0;
    % H
   % polyH = w4*eeH+w3;
   % polyH = eeH*polyH+w2;
   % polyH = eeH*polyH+w1;
   % polyH = eeH*polyH+w0;
    
    back_error_Pc(i) = rc/(polyc*norm(Xc(:,i)));
    %back_error_Pk(i)=rk/(polyk*norm(Xk(:,i)));
    %back_error_P1(i)=r1/(poly1*norm(X1(:,i)));
    %back_error_PH(i)=rH/(polyH*norm(XH(:,i)));
    %back_error_PG(i)=rG/(polyG*norm(XG(:,i)));
    
    vector_norm_ratioc(i) = norm(Vc(:,i))/norm(Xc(:,i));
    %vector_norm_ratiok(i) = norm(Vk(:,i))/norm(Xk(:,i));
    %vector_norm_ratio1(i) = norm(V1(:,i))/norm(X1(:,i));
end

%% EXACT COMPUTATIONS
if strcmp(symbolic,'yes')==1
    disp('Computing exact eigenpairs of the companion form')
    digits(64)
    A = vpa(sym(A,'f'));
    B = vpa(sym(B,'f'));
    %Compute eigenvalues, and right and left eigenvectors:
    [V,D] = eig(-B\A); %V contains the right eigenvectors
    [W,D2] = eig(-A'/B'); %W contains the left eigenvectors
    e = diag(D);
    e2 = diag(D2);
    %sort eigenvalues and eigenvectors according to absolute value:
    [~,ind] = sort(abs(e),'descend');
    e = e(ind);
    V = V(:,ind);
    [~,ind] = sort(abs(e2),'descend');
    %e2 = e2(ind);
    W = W(:,ind);

    disp('Computing exact eigenvectors of the matrix polynomial')
    X = sym(zeros(n,4*n)); %right eigenvectors
    Y = sym(zeros(n,4*n)); %left eigenvectors
    for i = 1:4*n
        Y(:,i)=W(1:n,i); %recover from first block
        if abs(e(i))>1
            X(:,i)=V(1:n,i); %recover from first block
        else
            X(:,i)=V(3*n+1:4*n,i); %recover from last block
        end
    end
% 
%     %% CONDITIONING OF THE EIGENVALUES OF THE POLYNOMIAL
%     disp('Conditioning of P')
%     condP = zeros(4*n,1);
%     for i=1:4*n
%         %Derivative
%         der = 4*e(i)*P4+3*P3;
%         der = der*e(i)+2*P2;
%         der = der*e(i)+P1;
%         if strcmp(wstring,'coeff')==1
%             w0 = n0;
%             w1 = n1;
%             w2 = n2;
%             w3 = n3;
%             w4 = n4;
%         else
%             w0 = nmax;
%             w1 = nmax;
%             w2 = nmax;
%             w3 = nmax;
%             w4 = nmax;
%         end
%        %Numerator
%        num = abs(e(i))*w4+w3;
%        num = abs(e(i))*num+w2;
%        num = abs(e(i))*num+w1;
%        num = abs(e(i))*num+w0;
% 
%        %Eigenvalue conditioning
%        condP(i)= num*norm(X(:,i))*norm(Y(:,i))/(abs(e(i))*abs(Y(:,i)'*der*X(:,i)));
%     end

    %% CONDITIONING OF THE EIGENVALUES OF D1 and DK 
    disp('conditioning of D1 and Dk')
    condDk = zeros(4*n,1);
    condD1 = zeros(4*n,1);
    R=zeros(4*n, 4*n);
    L=zeros(4*n,4*n);
    for i=1:4*n
        %right eigenvector formula
        R(:,i) = [e(i)^3*X(:,i) ; e(i)^2*X(:,i) ; e(i)*X(:,i);X(:,i)];
        %left eigenvector formula
        L(i,:) = [e(i)^3*Y(:,i)', e(i)^2*Y(:,i)', e(i)*Y(:,i)',Y(:,i)'];
        %Conditioning (normwise or coefficientwise)
        if strcmp(wstring,'coeff')==1
            numk = norm(R(:,i))*norm(L(i,:))*(nBk*abs(e(i))+nAk);
            num1 = norm(R(:,i))*norm(L(i,:))*(nB1*abs(e(i))+nA1);
        else
            numk = norm(R(:,i))*norm(L(i,:))*max([nBk nAk])*(1+abs(e(i)));
            num1 = norm(R(:,i))*norm(L(i,:))*max([nB1 nA1])*(1+abs(e(i)));
        end
        denk = abs(e(i))*abs(L(i,:)*Bk*R(:,i));
        den1 = abs(e(i))*abs(L(i,:)*B1*R(:,i));
        condDk(i) = numk/denk;
        condD1(i) = num1/den1;
    end
end

% % CONDITIONING OF THE EIGENVECTORS OF DK 
% disp('eigenvector conditioning of Dk')
% hk=zeros(4*n,1);
% %lowerDk=zeros(4*n,1);
% upperDk = zeros(4*n,1);
% 
% for i=1:4*n
%     hk(i)= 1- (norm(L(i,:)*Bk)*norm(R(:,i)))/abs(L(i,:)*Bk*R(:,i));
%     for k=1:4*n
%         if k~=i
%             %right eigenvector formula
%             %R = [e(k)^3*X(:,k) ; e(k)^2*X(:,k) ; e(k)*X(:,k);X(:,k)];
%             %left eigenvector formula
%             %L = [e(k)^3*Y(:,k)', e(k)^2*Y(:,k)', e(k)*Y(:,k)',Y(:,k)'];
%             s = norm(L(k,:))*norm(R(:,k))*(nAk+nBk*abs(e(i)))/(abs(e(i)-e(k))*abs(L(k,:)*Bk*R(:,k)));
%             upperDk(i)=(upperDk(i)+s);
%             %lowerDk(i)=upperDk(i)/(hk(i));
%         end
%     end
% end
% upperDk=4^(3/2)* eps*upperDk;
% %lowerDk=eps*lowerDk;

%% PLOTS

% EIGENVALUES
% semilogy(abs(e1),'o')

% CONDITION NUMBERS
%  semilogy(condP,'bo')
%  hold on
%  semilogy(condDk,'mx')
%  semilogy(condD1,'r*')

% BACKWARD ERRORS
% semilogy(back_error_Pk,'rx')
% hold on
% semilogy(back_error_P1,'b+')


% VECTOR NORM RATIOS
% semilogy(vector_norm_ratiok,'rx')
% hold on
% semilogy(vector_norm_ratio1,'b+')

% BACKWARD ERROR RATIOS
figure
semilogy(back_error_Pc./back_error_C,'rx')
%semilogy(back_error_Pk./back_error_Dk,'rx')
hold on
%semilogy(back_error_P1./back_error_D1,'b+')
%hold on
%semilogy(lowerDk, 'go')
%hold on
%hold on 
semilogy(vector_norm_ratioc,'bo')
hold on
%semilogy(vector_norm_ratiok, 'bo')
%semilogy(upperDk, 'yo')
legend('Pc/C')
%legend('Pk/Dk', 'P1/D1')
title('Degree 4 - backward')



% BACKWARD ERROR RATIOS - OTHER LINEARIZATIONS
%semilogy(back_error_PH./back_error_DH,'mx')
%semilogy(back_error_PG./back_error_DG,'g+')

%M4 = [P4 o o; P3 P4 o; P2 P3 P4];
%
