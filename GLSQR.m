function [X, REZ, ABS_E] = GLSQR(A, B, X_EXACT, maxit, X0, tol_stop)
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

R0=B-A*X0;

L(1,1)=norm(R0,'fro');
W(1:n,1:p) = R0/L(1,1);

L(1,2)=norm((A')*W(1:n,1:p),'fro');
S(1:n,1:p)=((A')*W(1:n,1:p))/L(1,2);


for k=1:maxit
    P=A* S(:,(p*(k-1)+1):(p*k))-L(k,k+1)*W(:,(p*(k-1)+1):(p*k));
    L(k+1,k+1)=norm(P,'fro');
    W(:,(p*k+1):(p*(k+1)))=P/L(k+1,k+1); 
    Q=(A')*W(:,(p*k+1):(p*(k+1)))-L(k+1,k+1)*S(:,(p*(k-1)+1):(p*k));    
    L(k+1,k+2)=norm(Q,'fro');
    S(:,(p*k+1):(p*(k+1)))=Q/L(k+1,k+2); 
    
    %%%
    T=eye(k+1);
    e1=T(:,1);
    
    y= L(1:k+1,2:(k+1)) \ ( L(1,1) * e1 );
    
    Xk=zeros(n,p);
    for i=1:k
        Xk=Xk+S(:,p*(i-1)+1:p*i)*y(i);    
    end
    
    REZ(k)=norm(B-A*Xk);    
    ABS_E(k)=norm(X_EXACT-Xk);
    
    if REZ(k)<tol_stop
       break 
    end    
    
end

X=Xk;

end
