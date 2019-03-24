%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Jednotlive bloky spustat s F9%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shaw, 2 prave strany
clear 
clc

tol_stop=1e-8;
tol_def=1e-12;
maxit=200;

n=100;

X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)]=shaw(n);
x(1:n,2)=sin(x(1:n,1));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%Blur, 2 prave strany
clear 
clc

tol_stop=1e-8;
tol_def=1e-12;
maxit=200;

n=16;
X0=zeros(n*n,2);
[A,B(1:n*n,1),x(1:n*n,1)] = blur (n);
x(1:n*n,2)=abs(sin(x(1:n*n,1)))*5;
B(1:n*n,2)=A*x(1:n*n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%Tomo, 2 prave strany
clear 
clc

tol_stop=1e-8;
tol_def=1e-12;
maxit=200;

n=10;
X0=zeros(n*n,2);
[A,B(1:n*n,1),x(1:n*n,1)] = tomo (n);
[A,B(1:n*n,2),x(1:n*n,2)] = tomo (n);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%deriv2, 2 prave strany, parameter 1,2,3
clear 
clc

tol_stop=1e-8;
tol_def=1e-12;
maxit=200;

n=32;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = deriv2 (n,3);
%B(1:n,1) = B(1:n,1) + 1e-3*randn(size(B(1:n,1))); %%%
x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%deriv2, 2 prave strany
clear 
clc

tol_stop=1e-10;
tol_def=1e-12;
maxit=200;

n=32;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = foxgood (n);

x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%gravity, 2 prave strany
clear 
clc

tol_stop=1e-12;
tol_def=1e-12;
maxit=200;

n=100;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = gravity (n,1);

x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%heat, 2 prave strany, kappa=5 well-possed, kappa=1 ill-possed
%velmi zaujimave
clear 
clc

tol_stop=1e-12;
tol_def=1e-12;
maxit=200;

n=100;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = heat (n,1);

x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%phillips , 2 prave strany
clear 
clc

tol_stop=1e-10;
tol_def=1e-12;
maxit=200;

n=100;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = phillips(n);

x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

%wing , 2 prave strany --nespojite riesenie
clear 
clc

tol_stop=1e-12;
tol_def=1e-12;
maxit=200;

n=100;
X0=zeros(n,2);

[A,B(1:n,1),x(1:n,1)] = wing(n);

x(1:n,2) = x(1:n,1) + 1e-3*randn(size(x(1:n,1)));
B(1:n,2)=A*x(1:n,2);

X_EXACT=x;

testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)
