function testD(A, B, X_EXACT, maxit, X0, tol_stop, tol_def)
%%%

disp('BGMRES:')
tic
[X1, REZ1, ABS_E1] = BGMRESdef(A, B, X_EXACT, maxit, X0, tol_stop, tol_def);
toc

subplot(2,4,1), 
semilogy(REZ1), title('residuum BGMRES'),

subplot(2,4,5), 
semilogy(ABS_E1), title('abs. error BGMRES'),
%%%

disp('GGMRES:')
tic
[X2,REZ2]=GGMRES(A,B,tol_stop,maxit,X0);
toc

subplot(2,4,2)
semilogy(REZ2), title('residuum GGMRES'),

%%%

disp('BLSQR:')
tic
[X3, REZ3, ABS_E3] = BLSQRdef(A, B, X_EXACT, maxit, X0, tol_stop, tol_def);
toc

subplot(2,4,3)
semilogy(REZ3), title('residuum BLSQR'),

subplot(2,4,7)
semilogy(ABS_E3), title('abs. error BLSQR'),

%%%

disp('GLSQR:')
tic
[X4,REZ4]=GLSQR(A,B,tol_stop,maxit,X0);
toc

subplot(2,4,4)
semilogy(REZ4), title('residuum GLSQR'),

end