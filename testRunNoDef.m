function testRunNoDef(A, B, X_EXACT, maxit, X0, tol_stop)
%%%

disp('BGMRES:')
tic
[X1, REZ1, ABS_E1] = BGMRES(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,1), 
semilogy(REZ1), title('residuum BGMRES'),

subplot(2,6,7), 
semilogy(ABS_E1), title('abs. error BGMRES'),
%%%

disp('GGMRES:')
tic
[X2,REZ2, ABS_E2]=GGMRES(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,2)
semilogy(REZ2), title('residuum GGMRES'),

subplot(2,6,8), 
semilogy(ABS_E2), title('abs. error GGMRES'),

%%%

disp('BLSQR:')
tic
[X3, REZ3, ABS_E3] = BLSQR(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,3)
semilogy(REZ3), title('residuum BLSQR'),

subplot(2,6,9)
semilogy(ABS_E3), title('abs. error BLSQR'),

%%%

disp('GLSQR:')
tic
[X4,REZ4, ABS_E4]=GLSQR(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,4)
semilogy(REZ4), title('residuum GLSQR'),

subplot(2,6,10)
semilogy(ABS_E4), title('abs. error GLSQR'),

%%%

disp('CGMRES:')
tic
[X6,REZ5, ABS_E5]=CGMRES(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,5)
semilogy(REZ5), title('residuum CGMRES'),

subplot(2,6,11)
semilogy(ABS_E5,'*'), title('abs. error CGMRES'),

%%%

disp('CLSQR:')
tic
[X6,REZ6, ABS_E6]=CLSQR(A, B, X_EXACT, maxit, X0, tol_stop);
toc

subplot(2,6,6)
semilogy(REZ6), title('residuum CLSQR'),

subplot(2,6,12)
semilogy(ABS_E6,'*'), title('abs. error CLSQR'),


end