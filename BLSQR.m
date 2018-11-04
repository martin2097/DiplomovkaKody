function [X,REZ] = BLSQR(A,B,tol,maxit,X0)
n = size(A,1);
p = size(B,2);

R0=B-A*X0;
[S1,D1] = qr(R0,0);
S(:,1:p) = S1;
Wr=A*S(:,1:p);
[W(:,1:p),L(1:p,1:p)]=qr(Wr,0);

for k=1:maxit
    Fr=(A')* W(:,(p*(k-1)+1):(p*k))-S(:,(p*(k-1)+1):(p*k))*((L(p*(k-1)+1:p*k,p*(k-1)+1:p*k))');
    Fr=Fr-S(:,1:p*k)*((S(:,1:p*k)')*Fr);
    if k<maxit
       [S(:,(p*k+1):(p*(k+1))),Rj]=qr(Fr,0);
       
%norm(S(:,(p*k+1):(p*(k+1)))'*S(:,(p*k+1):(p*(k+1)))-eye(p,p)) %test 1)
%norm(S'*S-eye((k+1)*p,(k+1)*p)) %test 2)



       
       L((p*(k-1)+1):(p*k),(p*k+1):(p*(k+1)))=Rj';
       Wr=A*S(:,(p*k+1):(p*(k+1)))-W(:,(p*(k-1)+1):(p*k))*L((p*(k-1)+1):(p*k),(p*k+1):(p*(k+1)));
       Wr=Wr-W(:,1:p*k)*((W(:,1:p*k)')*Wr);
       [W(:,(p*k+1):(p*(k+1))),L(p*k+1:p*(k+1),p*k+1:p*(k+1))]=qr(Wr,0);
       
%norm(W(:,(p*k+1):(p*(k+1)))'*W(:,(p*k+1):(p*(k+1)))-eye(p,p)) %test 1)
%norm(W'*W-eye((k+1)*p,(k+1)*p))  %test 2)


    
    %%%
    
    E1=eye((k+1).*p,p);
    Y=L(1:(k+1)*p,p+1:(k+1)*p)\(E1*L(1:p,1:p));
    Xk=W(:,1:p*k)*Y;
    REZ(k)=norm(B-A*Xk);    
   

    if (REZ(k)<tol)
       break 
    end
    end
    
    
end

   X=Xk;
   % norm(A'*W(:,1:(maxit-1)*p)-S*L(1:(maxit-1)*p,1:maxit*p)') %test 3


end