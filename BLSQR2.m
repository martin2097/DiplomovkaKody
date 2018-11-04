function [X,REZ] = BLSQR2(A,B,X0,m)
n = size(A,1);
p = size(B,2);

W(:,1:p)=zeros(n,p); %W0
R0=B-A*X0;
[S1,D1] = qr(R0,0);

S(:,1:p) = S1;
L(1:p,1:p)=D1;

for k=1:m
    P=A* S(:,(p*(k-1)+1):(p*k))-W(:,(p*(k-1)+1):(p*k))*L((p*(k-1)+1):(p*k),(p*(k-1)+1):(p*k));
    [WK,GK]=qr(P',0);
    W(:,(p*k+1):(p*(k+1)))=WK'; %POZOR! Do W davam aj blok W0 teda Wm=[W0,...,Wm]
    L((p*(k-1)+1):(p*k),(p*k+1):(p*(k+1)))=GK';
    Q=A*W(:,(p*k+1):(p*(k+1)))-S(:,(p*(k-1)+1):(p*k))*L((p*(k-1)+1):(p*k),(p*k+1):(p*(k+1)));    
    [SK,DK]=qr(Q,0);    
    S(:,(p*k+1):(p*(k+1)))=SK; %Tu je Sm=[S1,...,S(m+1)]
    L((p*k+1):(p*(k+1)),(p*k+1):(p*(k+1)))=DK; %Rovnako matica Lk je horna trojuholnikova
    
    
end












end