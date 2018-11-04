function [X,REZ] = GGMRES(A,B,tol,maxit,X0)
n = size(A,1);
p = size(B,2);


R0=B-A*X0;
V1=R0/norm(R0,'fro');
%%%

%H = zeros(m);
%U = zeros(n,m*p);

U(1:n,1:p) = V1;
for j=1:maxit
    for i=1:j
    H(i,j)=trace(U(1:n,p*(i-1)+1:p*i)'*A*U(1:n,p*(j-1)+1:p*j));
    end
    S=zeros(n,p);
    for i=1:j
    S=S+U(1:n,p*(i-1)+1:p*i)*H(i,j);    
    end    
    W=A*U(1:n,p*(j-1)+1:p*j)-S;
    H(j+1,j)=sqrt(trace(W'*W));
    U(1:n,p*(j)+1:p*(j+1))=1/H(j+1,j)*W;
    
    %%%
    T=eye(j+1);
    e1=T(:,1);

    y= H \ ( norm(R0,'fro') * e1 );
    

    %%%vysledok%%%
    Xj=X0;
    for i=1:j
        Xj=Xj+U(:,p.*(i-1)+1:p.*i).*y(i);    
    end
    
    REZ(j)=norm(B-A*Xj);
    
    if REZ(j)<tol
       break 
    end    
end

X=Xj;
%%%




end