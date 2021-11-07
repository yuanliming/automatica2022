syms ki1 kp1 kd1 ki2 kp2 kd2 twokm_L
K=[ki1   0       kp1        0  kd1    0;
   0   ki2         0      kp2    0  kd2;
   0     0  -twokm_L  twokm_L    0    0;]; % structured feedback matrix K of example 5.1
K1=[ki1 kp1 kd1];       
K2=[ki2 kp2 kd2];       
K3=twokm_L;             
[r,TC,TR]=algo1_facK(K,{K1,K2,K3}); % factorize K=TCKDTR, where KD is block diagonal matrix with diagonal blocks of K1, K2, K3; TC, TR are constant matrices.
KD=blkdiag(kron(eye(r(1)),K1),kron(eye(r(2)),K2),kron(eye(r(3)),K3));
disp(KD)
disp(TC)
disp(TR)
