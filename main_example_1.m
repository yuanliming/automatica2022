syms k1 k2 k3 
K=[k1    -k1   0;
   -k1 k1+k2  k3;]; 
K1=k1;      
K2=[k2 k3];       
[r,TC,TR]=algo1_facK(K,{K1,K2}); 
KD=blkdiag(kron(eye(r(1)),K1),kron(eye(r(2)),K2));
disp(KD)
disp(TC)
disp(TR)
