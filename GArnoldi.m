function [U,H] = GArnoldi(A,B,U1,m)
%Zovöeobecnen˝ Arnoldiho proces podæa Algoritmu 2.1 z Ël·nku
%Vstup:     A - matica n x n
%           m - poËet iter·ciÌ
%           p - n·sobnosù pravej strany (öÌrka bloku)
%V˝stup:    U - blokovo unit·rna matica
%           H - blokovo Hessenbergova matica
n = size(A,1);
p = size(B,2);
H = zeros(m);
U = zeros(n,m*p);

U(:,1:p) = U1;
for j=1:m
    for i=1:j
    H(i,j)=trace(U(:,p*(i-1)+1:p*i)'*A*U(:,p*(j-1)+1:p*j));
    end
    S=zeros(n,p);
    for i=1:j
    S=S+U(:,p*(i-1)+1:p*i)*H(i,j);    
    end    
    W=A*U(:,p*(j-1)+1:p*j)-S;
    H(j+1,j)=sqrt(trace(W'*W));
    U(:,p*(j)+1:p*(j+1))=1/H(j+1,j)*W;
       
end

end