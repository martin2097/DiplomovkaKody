function [U,H] = BArnoldi(A,B,U1,m)
%Blokový Arnoldiho proces pod¾a Algoritmu 6.21 zo Saada
%Vstup:     A - matica n x n
%           m - poèet iterácií
%           B - štartovací blok
%Výstup:    U - blokovo unitárna matica
%           H - blokovo Hessenbergova matica
n = size(A,1);
p = size(B,2);
H = zeros(m*p);
U = zeros(n,m*p);

U(:,1:p) = U1;
for j=1:m
    for i=1:j
    H((p*(i-1)+1):(i*p),(p*(j-1)+1):(j*p))=U(:,p*(i-1)+1:p*i)'*A*U(:,p*(j-1)+1:p*j);
    end
    S=zeros(n,p);
    for i=1:j
    S=S+U(:,p*(i-1)+1:p*i)*H((p*(i-1)+1):(i*p),(p*(j-1)+1):(j*p));    
    end    
    W=A*U(:,p*(j-1)+1:p*j)-S;
    [U(:,p*(j)+1:p*(j+1)),H((p*((j+1)-1)+1):((j+1)*p),(p*(j-1)+1):(j*p))] = qr(W,0);
       
end




end