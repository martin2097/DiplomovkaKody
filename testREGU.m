%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Jednotlive bloky spustat s F9%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shaw, 2 prave strany
clear 
clc

tol=1e-8;
maxit=200;

n=100;

X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)]=shaw(n);
x(1:n,2)=sin(x(1:n,1));
B(1:n,2)=A*x(1:n,2);

test(A,B,tol,maxit,X0)

%Blur, 2 prave strany
clear 
clc

tol=1e-8;
maxit=200;

n=16;
X0=zeros(n*n,2);
[A,B(1:n*n,1),x(1:n*n,1)] = blur (n);
x(1:n*n,2)=abs(sin(x(1:n*n,1)))*5;
B(1:n*n,2)=A*x(1:n*n,2);

test(A,B,tol,maxit,X0)

%Tomo, 2 prave strany
clear 
clc

tol=1e-6;
maxit=200;

n=10;
X0=zeros(n*n,2);
[A,B(1:n*n,1),x(1:n*n,1)] = tomo (n);
[A,B(1:n*n,2),x(1:n*n,2)] = tomo (n);

test(A,B,tol,maxit,X0)
