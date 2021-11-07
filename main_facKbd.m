syms k1 k2 k3 k4 k5 k5 k6
K=[k1 k2 k3 k1 k2 k3;
   k4 k5 k6 k4 k5 k6;
   0  0  0  k1 k2 k3;
   0  0  0  k4 k5 k6]; % structured feedback matrix K of example 5.3
K1=[k1 k2 k3;
    k4 k5 k6]; % the desired block of control parameters.
[r,TC,TR]=algo1_facK(K,{K1}); 
KD=blkdiag(kron(eye(r(1)),K1));
disp(KD)
disp(TC)
disp(TR)