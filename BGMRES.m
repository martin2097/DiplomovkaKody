function [X,REZ] = BGMRES(A,B,tol,maxit,X0)
n = size(A,1);
p = size(B,2);

R0=B-A*X0;
[V1,R] = qr(R0,0);

U(:,1:p) = V1;
for j=1:maxit
    for i=1:j
    H((p*(i-1)+1):(i*p),(p*(j-1)+1):(j*p))=U(:,p*(i-1)+1:p*i)'*A*U(:,p*(j-1)+1:p*j);
    end
    S=zeros(n,p);
    for i=1:j
    S=S+U(:,p*(i-1)+1:p*i)*H((p*(i-1)+1):(i*p),(p*(j-1)+1):(j*p));    
    end    
    W=A*U(:,p*(j-1)+1:p*j)-S;
    [U(:,p*(j)+1:p*(j+1)),H((p*((j+1)-1)+1):((j+1)*p),(p*(j-1)+1):(j*p))] = qr(W,0);
    %%%
    E1=eye((j+1).*p,p);
    Y=H\(E1*R);
    Xk=X0+U(:,1:(j.*p))*Y;
               
    REZ(j)=norm(B-A*Xk);
    if REZ(j)<tol
       break 
    end
end

X=Xk;

end