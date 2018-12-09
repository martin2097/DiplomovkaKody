function test(A,B,tol,maxit,X0)
%%%

disp('BGMRES:')
tic
[X1,REZ1]=BGMRES(A,B,tol,maxit,X0);
toc

subplot(2,3,1), 
semilogy(REZ1), title('residuum BGMRES'),

%%%

disp('GGMRES:')
tic
[X2,REZ2]=GGMRES(A,B,tol,maxit,X0);
toc

subplot(2,3,2)
semilogy(REZ2), title('residuum GGMRES'),

%%%

disp('CGMRES:')
tic
[X3,REZ3] = CGMRES(A,B,tol,maxit);
toc

subplot(2,3,3)
semilogy(REZ3), title('residuum gmres buid-in'),

%%%

disp('BLSQR:')
tic
[X4,REZ4]=BLSQR(A,B,tol,maxit,X0);
toc

subplot(2,3,4)
semilogy(REZ4), title('residuum BLSQR'),

%%%

disp('GLSQR:')
tic
[X5,REZ5]=GLSQR(A,B,tol,maxit,X0);
toc

subplot(2,3,5)
semilogy(REZ5), title('residuum GLSQR'),

%%%

disp('CLSQR:')
tic
[X6,REZ6] = CLSQR(A,B,tol,maxit);
toc

subplot(2,3,6)
semilogy(REZ6), title('residuum lsqr buid-in'),

end