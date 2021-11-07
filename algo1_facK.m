function [r,TC,TR] = algo1_facK(K,k)
%% factorize K=TCKDTR, 
%% KD is block diagonal matrix with diagonal blocks of k={K1,K2,..,KN}
%% TC, TR are constant matrices. see Theorem 1.
%% any question/adivice contact yuanliming@nimte.ac.cn, chensilu@nimte.ac.cn
r = []; TC = [];TR = [];
[m,n] = size(K);
for i = 1:length(k)
bar_Ki = pt(K,k{i}); % eq.(9)
Ti = obT(bar_Ki,k{i}); % eq.(10)
[mi,ni] = size(k{i});
til_Ti = rar(Ti,m,mi); % eq.(17)
ri = rank(til_Ti);
[Fi,Gi] = fr(til_Ti);
TRi = reshape(Fi,n,ri*ni)'; % eq.(19)
TCi = reshape(Gi',m,ri*mi); % eq.(20)
%KD=kron(eye(r),k);
r = [r ri];
TC = [TC TCi];
TR = [TR;TRi];
end
end
function [barKi] = pt(K,Ki) 
%% partition a matrix (Definition 3)
[m,n] = size(K);
[mi,ni] = size(Ki);
barKi = zeros(m,n);
vecKi = reshape(Ki,mi*ni,1);
for i = 1:length(vecKi)
subKi = K-subs(K,vecKi(i),0);
barKi = barKi+subKi;
end
end
function [Ti] = obT(bar_Ki,K_i)
%% obtain Ti by sovling eq.(18)
[m,n] = size(bar_Ki);
[mi,ni] = size(K_i);
vecKi = reshape(K_i,mi*ni,1); % vectorize Ki
vecbar_Ki = reshape(bar_Ki,m*n,1); % vectorize \bar Ki
Ti = [];
I = eye(length(vecKi));
for i = 1:length(vecKi) % obtain Ti
Tt = subs(vecbar_Ki,vecKi,I(:,i));
Ti = [Ti Tt];
end
end
function [til_Ti] = rar(Ti,m,mi) 
%% rearrange a matrix (Definition 5)
til_Ti = [];
[mT,nT] = size(Ti);
bk_T = mat2cell(Ti,repmat(m,1,mT/m),repmat(mi,1,nT/mi)); % partition T into a block matrix with each block (sub-matrice) of size [m mi].
nbk = mT/m*nT/mi; % there are nbk sub-matrices in total
for i = 1:nbk
vecTi = reshape(bk_T{i},m*mi,1); % vectorize every sub-matrice
til_Ti = [til_Ti;vecTi']; % stacking the vectorized sub-matrices from top to bottom
end
til_Ti = double(til_Ti);
end
function [B,C] = fr(A) 
%% rank facotriztion of matrix A=BC (Definition 6)
[E,p] = rref(A);
B = A(:,p);
C = E(1:length(p),:);
end