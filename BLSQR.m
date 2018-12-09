function [X,REZ] = BLSQR(A,B,tol,maxit,X0)
n = size(A,1);
p = size(B,2);
sw=0;

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
    [V(1:n,p*k+1:p*(k+1)),T(p*k+1:p*(k+1),p*(k+1)+1:p*(k+2))]=qr(Zk,0);
       
    %%%
    
    E1=eye((k+1).*p,p);
    Y=T(1:(k+1)*p,p+1:(k+1)*p)\(E1*T(1:p,1:p));
    Xk=V(:,1:p*k)*Y;
    REZ(k)=norm(B-A*Xk);    
    
   
    if rank(V,10e-14)<(k+1)*p 
       sw=1; 
    end
    
    if (REZ(k)<tol)
       break 
    end
    end
    
    X=Xk; 
    
    %V(:,1:p*k)'*V(:,1:p*k)
    no=norm(V(:,1:p*k)'*V(:,1:p*k)-eye(size(V(:,1:p*k),2)));
    suc=0;
    for i=1:k
        suc=suc+norm(V(:,p*(i-1)+1:p*i)'*V(:,p*(i-1)+1:p*i)-eye(p));
    end
    
    disp(['Strata ortogonality celej matice je ',num2str(no),' a po blokoch ',num2str(suc)]);
    
    if sw==1 
       disp('Deflace!')
    end
    
end

  
   
