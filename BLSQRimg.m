function [X, REZ, ABS_E] = BLSQRimg(B, X_EXACT, maxit, X0, tol_stop, tol_def)
n = size(B,1);
p = size(B,2);
sigma = 1;

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

Usize = zeros(1,maxit+2);
Vsize = zeros(1,maxit+2);

R0=B-imgaussfilt(X0,sigma);
[U1,R1] = qr(R0,0);

%test deflacie v U1
it=0;
for i=1:p
    it=it+1;
    if norm(U1(:,it))<tol_def
        U1(:,it)=[];
        R1(it,:)=[];
        it=it-1;
    end
end
Usize(1,2)= size(U1,2);
U(:,1:Usize(1,2)) = U1;
T(1:Usize(1,2),1:p)=R1;

Z1=(A')*U(:,1:Usize(1,2));

[V1,L1]=qr(Z1,0);

%test deflacie vo V1
it=0;
for i=1:Usize(1,2)
    it=it+1;
    if norm(V1(:,it))<tol_def
        V1(:,it)=[];
        L1(it,:)=[];
        it=it-1;
    end
end
Vsize(1,2)= p;
Vsize(1,3)= p + size(V1,2);

V(:,1:Vsize(1,3)-Vsize(1,2))=V1;
T(1:Usize(1,2),p+1:Vsize(1,3))=L1';

for k=1:maxit
    Wk=A*V(:,Vsize(1,k+1)+1-Vsize(1,2):Vsize(1,k+2)-Vsize(1,2))-U(:,Usize(1,k)+1:Usize(1,k+1))*T(Usize(1,k)+1:Usize(1,k+1),Vsize(1,k+1)+1:Vsize(1,k+2));
    [U_k,R_k]=qr(Wk,0);
    
    %test deflacie v U_k
    it=0;
    for i=1:(Vsize(1,k+2)-Vsize(1,k+1))
        it=it+1;
        if norm(U_k(:,it))<tol_def
            U_k(:,it)=[];
            R_k(it,:)=[];
            it=it-1;
        end
    end
    Usize(1,k+2)=Usize(1,k+1)+size(U_k,2);
    U(:,Usize(1,k+1)+1:Usize(1,k+2)) = U_k;
    
    T(Usize(1,k+1)+1:Usize(1,k+2),Vsize(1,k+1)+1:Vsize(1,k+2))=R_k;
    
    Zk=A'*U(:,Usize(1,k+1)+1:Usize(1,k+2))-V(:,Vsize(1,k+1)+1-Vsize(1,2):Vsize(1,k+2)-Vsize(1,2))*T(Usize(1,k+1)+1:Usize(1,k+2),Vsize(1,k+1)+1:Vsize(1,k+2))';
    [V_k,L_k]=qr(Zk,0);
    
    %test deflacie v V_k
    it=0;
    for i=1:(Usize(1,k+2)-Usize(1,k+1))
        it=it+1;
        if norm(V_k(:,it))<tol_def
            V_k(:,it)=[];
            L_k(it,:)=[];
            it=it-1;
        end
    end
    
    Vsize(1,k+3)= Vsize(1,k+2) + size(V_k,2);

    V(:,Vsize(1,k+2)+1-Vsize(1,2):Vsize(1,k+3)-Vsize(1,2))=V_k;
    T(Usize(1,k+1)+1:Usize(1,k+2),Vsize(1,k+2)+1:Vsize(1,k+3))=L_k';
    
    %%%
    
    E1=eye((k+1).*p,p);
    Y=T(1:Usize(k+2),Vsize(2)+1:Vsize(k+2))\(E1*T(1:Usize(1,2),1:p));
    
    Xk=V(:,1:Vsize(1,k+2)-Vsize(1,2))*Y;
    REZ(k)=norm(B-A*Xk);    
    ABS_E(k)=norm(X_EXACT-Xk);
   
    if (REZ(k)<tol_stop) || (Usize(1,k+2) - Usize(1,k+1) == 0) || (Vsize(1,k+3) - Vsize(1,k+2) == 0)
       break 
    end
end
   
X=Xk; 

end

  
   
