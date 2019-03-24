function [X, REZ, ABS_E] = BLSQR(A, B, X_EXACT, maxit, X0, tol_stop)
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
[U1,R1] = qr(R0,0);
T(1:p,1:p)=R1;
U(:,1:p) = U1;
Z1=(A')*U(:,1:p);
[V(:,1:p),L1]=qr(Z1,0);
T(1:p,p+1:2*p)=L1';

for k=1:maxit
    Wk=A*V(:,p*(k-1)+1:p*k)-U(:,p*(k-1)+1:p*k)*T(p*(k-1)+1:p*k,p*k+1:p*(k+1));
    [U(1:n,p*k+1:p*(k+1)),T(p*k+1:p*(k+1),p*k+1:p*(k+1))]=qr(Wk,0);
    
    Zk=A'*U(1:n,p*k+1:p*(k+1))-V(:,p*(k-1)+1:p*k)*T(p*k+1:p*(k+1),p*k+1:p*(k+1))';
    [V(1:n,p*k+1:p*(k+1)),L_k]=qr(Zk,0);
    T(p*k+1:p*(k+1),p*(k+1)+1:p*(k+2))=L_k';
    %%%
    
    E1=eye((k+1).*p,p);
    Y=T(1:(k+1)*p,p+1:(k+1)*p)\(E1*T(1:p,1:p));
    Xk=V(:,1:p*k)*Y;
    
    REZ(k)=norm(B-A*Xk);    
    ABS_E(k)=norm(X_EXACT-Xk);
       
    if (REZ(k)<tol_stop)
       break 
    end
end
  
X=Xk;   
    
end

  
   
