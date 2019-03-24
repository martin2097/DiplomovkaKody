%cluster 2,-4
clear 
clc

tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=randn(n,p);

A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 10*randn(n,n)/n;      % cluster 2,-4

B=A*X_EXACT;

%plot(eig(A),'*')                   

testRunNoDef(A, B, X_EXACT, maxit, X0, tol_stop)

%cluster 2, sirsie okolie
clear 
clc

tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);

A=2*eye(n) + 10*randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci

%plot(eig(A),'*')                   

testRunNoDef(A, B, X_EXACT, maxit, X0, tol_stop)

%cluster 2, uzsie okolie
clear 
clc

tol_stop=1e-10;
maxit=200;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=zeros(n,p);

B=randn(n,p);
A=2*eye(n) + randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci

%plot(eig(A),'*')                   

testRunNoDef(A, B, X_EXACT, maxit, X0, tol_stop)

%interesting
clear 
clc
tol_def=1e-4;
tol_stop=1e-8;
maxit=400;

n=300;
p=8;

X0=zeros(n,p);
X_EXACT=randn(n,p);


A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 40*randn(n,n)/n;    % matice s vl. c. 2 porusena o nahodnou perturbaci

B=A*X_EXACT;

%plot(eig(A),'*')                   

testRunNoDef(A, B, X_EXACT, maxit, X0, tol_stop)
