function [X, REZ, ABS_E] = CLSQR(A, B, X_EXACT, maxit, X0, tol_stop)
n = size(A,1);
p = size(B,2);

if ~exist('tol_stop','var')
     tol_stop = 1e-10;
 end
 
 if ~exist('X0','var')
     X0=zeros(n,p);
 end
 
 if ~exist('maxit','var')
     maxit = 200;
 end
 
 if ~exist('X_EXACT','var')
     X_EXACT=zeros(n,p);
 end

REZp=zeros(maxit,1);
iter=0;
maxitn=min(size(A,1),maxit);

    for k=1:p
    
    [Xk,flag,relres,iterk,REZk] = lsqr(A,B(:,k),tol_stop,maxitn);    
    X(1:n,k)=Xk;
    iter=max(iterk,iter);
    REZp(1:maxit)=REZp+[REZk;(ones(maxit-iterk-1,1)*REZk(iterk+1,1))];

    end

REZ=REZp(1:iter+1);
ABS_E=norm(X_EXACT-X);

end