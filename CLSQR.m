function [X,REZ] = CLSQR(A,B,tol,maxit)
n = size(A,1);
p = size(B,2);

REZp=zeros(maxit,1);
iter=0;
maxitn=min(size(A,1),maxit);

    for k=1:p
    
    [Xk,flag,relres,iterk,REZk] = lsqr(A,B(:,k),tol,maxitn);    
    X(1:n,k)=Xk;
    iter=max(iterk,iter);
    REZp(1:maxit)=REZp+[REZk;(ones(maxit-iterk-1,1)*REZk(iterk+1,1))];

    end

REZ=REZp(1:iter+1);




end