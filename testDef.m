%cluster 2,-4
clear 
clc
tol_def=1e-4;
tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);
A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 10*randn(n,n)/n;      % cluster 2,-4

%plot(eig(A),'*')                   

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%cluster 2, sirsie okolie
clear 
clc
tol_def=1e-4;
tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);

A=2*eye(n) + 10*randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci

%plot(eig(A),'*')                   

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%cluster 2, uzsie okolie
clear 
clc
tol_def=1e-4;
tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);
A=2*eye(n) + randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci

%plot(eig(A),'*')                   

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%interesting
clear 
clc
tol_def=1e-4;
tol_stop=1e-8;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);
A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 40*randn(n,n)/n;    % matice s vl. c. 2 porusena o nahodnou perturbaci

%plot(eig(A),'*')                   

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%test def
%cluster 2,-4
clear 
clc
tol_def=1e-4;
tol_stop=1e-10;
maxit=200;

n=300;
p=2;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B(1:n,1)=randn(n,1);


A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 10*randn(n,n)/n;      % cluster 2,-4
B(1:n,2)=A*B(1:n,1);
%plot(eig(A),'*')                   

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

