clear, clc

n=300;
tol=1e-8;
maxit=200;
p=8;

B=randn(n,p);
%B=[B,(B + randn(n,1)/n)];

%A=2*eye(n) + randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci
%A=2*eye(n) + 10*randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci
A=blkdiag(2*eye(n/2), -4*eye(n/2)) + 10*randn(n,n)/n;      % matice s vl. c. 2 porusena o nahodnou perturbaci

X0=zeros(n,p);
%plot(eig(A),'*')                   % -> A ma klastrovana vl. cisla



%%%

test(A,B,tol,maxit,X0)




