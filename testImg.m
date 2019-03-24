clear 
clc

A=im2double(imread('pumpkins.tif'));
B = imgaussfilt(A, 1);

X_EXACT=A;
maxit=10;

[X, REZ, ABS_E] = BGMRESimg(B, X_EXACT, maxit);

disp('BGMRESdef:')
subplot(1,3,1),
imshow(X),

subplot(1,3,2), 
semilogy(REZ), title('reziduum defBGMRES'),

subplot(1,3,3), 
semilogy(ABS_E), title('absolutna chyba defBGMRES'),

