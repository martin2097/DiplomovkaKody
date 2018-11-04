function [X,REZ] = CGMRES(A,B,tol,maxit)
n = size(A,1);
p = size(B,2);

REZp=zeros(maxit,1);
iter=0;
restart=min(size(A,1),maxit);

    for k=1:p
    
	[Xk,flag,relres,iterk,REZk] = gmres(A,B(:,k),restart,tol,restart);
    X(1:n,k)=Xk;
    iter=max(iterk(1,2),iter);
    REZp(1:maxit)=REZp+[REZk;(ones(maxit-iterk(1,2)-1,1)*REZk(iterk(1,2)+1,1))];

    end

REZ=REZp(1:iter+1);




end