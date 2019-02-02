clear 
clc

tol=1e-12;
maxit=200;

n=16;
X0=zeros(n*n,2);
[A,B(1:n*n,1),x(1:n*n,1)] = blur (n);
x(1:n*n,2)=abs(sin(x(1:n*n,1)))*5;
B(1:n*n,2)=A*x(1:n*n,2);


%%%%%%%%%%%%%%%%%%

disp('BGMRES:')
tic
[X1,REZ1]=BGMRES(A,B,tol,maxit,X0);
toc

subplot(2,3,1), 
semilogy(REZ1), title('reziduum BGMRES'),

%%%%%
disp('BGMRESdef:')
tic
[X2, REZ2, ABS_E2] = BGMRESdef(A, B,x, maxit, X0, tol);
toc

subplot(2,3,2), 
semilogy(REZ2), title('reziduum defBGMRES'),

subplot(2,3,3), 
semilogy(ABS_E2), title('absolutna chyba defBGMRES'),

%%%%
disp('BLSQR:')
tic
[X4,REZ4]=BLSQR(A,B,tol,maxit,X0);
toc

subplot(2,3,4)
semilogy(REZ4), title('reziduum BLSQR'),

%%%%%
disp('BLSQRdef:')
tic
[X5, REZ5, ABS_E5] = BLSQRdef(A, B, x, maxit, X0, tol);
toc

subplot(2,3,5), 
semilogy(REZ5), title('reziduum defBLSQR'),

subplot(2,3,6), 
semilogy(ABS_E5), title('absolutna chyba defBLSQR'),



