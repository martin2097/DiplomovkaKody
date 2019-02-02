function [X, REZ, ABS_E] = BGMRESdef(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)

n = size(A,1);
p = size(B,2);

 if ~exist('tol_def','var')
     tol_def = 1e-12;
 end
 
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
[V1,R] = qr(R0,0);

Usize = zeros(1,maxit+1);

%test deflacie v V1
it=0;
for i=1:p
    it=it+1;
    if norm(V1(:,it))<tol_def
        V1(:,it)=[];
        R(it,:)=[];
        it=it-1;
    end
end
Usize(1,2)= size(V1,2);
U(:,1:Usize(1,2)) = V1;

for j=1:maxit
    for i=1:j
    H((Usize(1,i)+1):(Usize(1,i+1)),(Usize(1,j)+1):(Usize(1,j+1)))=U(:,(Usize(1,i)+1):(Usize(1,i+1)))'*A*U(:,(Usize(1,j)+1):(Usize(1,j+1)));
    end
    S=zeros(n,Usize(1,i+1)-Usize(1,i));
    for i=1:j
    S=S+U(:,(Usize(1,i)+1):(Usize(1,i+1)))*H((Usize(1,i)+1):(Usize(1,i+1)),(Usize(1,j)+1):(Usize(1,j+1)));    
    end    
    W=A*U(:,(Usize(1,j)+1):(Usize(1,j+1)))-S;
    [Unew, Hnew] = qr(W,0);
        
    H((Usize(1,j+1)+1):(2*Usize(1,j+1))-Usize(1,j),(Usize(1,j)+1):(Usize(1,j+1))) = Hnew;
    
    %test deflacie vo V_k
    it=0;
    for i=1:(Usize(1,j+1)-Usize(1,j))
        it=it+1;
        if norm(Unew(:,it))<tol_def
           Unew(:,it)=[];
           H(Usize(1,j+1)+it,:) = [];
            it=it-1;
        end
    end
    Usize(1,j+2)= Usize(1,j+1)+size(Unew,2);
    U(:,Usize(1,j+1)+1:Usize(1,j+2)) = Unew;
    %%%
    E1=eye((j+1).*p,p);
    Y=H\(E1*R);
    Xk=X0+U(:,1:Usize(1,j+1))*Y;
         
    REZ(j)=norm(B-A*Xk);
    ABS_E(j)=norm(X_EXACT-Xk);
    
    if (REZ(j)<tol_stop) || (Usize(1,j+2) - Usize(1,j+1) == 0)
       break 
    end
end

X=Xk;

end