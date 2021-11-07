syms k11 k12 k23 k34 k45 k55
K=[k11+k12      -k12         0        0        0;
      -k12   k12+k23      -k23        0        0;
         0      -k23   k23+k34     -k34        0;
         0         0      -k34  k34+k45     -k45;
         0         0        0      -k45  k45+k55;]; % structured feedback matrix K of example 5.2 
K1=k11; K2=k12; K3=k23; K4=k34; K5=k45; K6=k55; 
[r,TC,TR]=algo1_facK(K,{K1,K2,K3,K4,K5,K6});
KD=blkdiag(kron(eye(r(1)),K1),kron(eye(r(2)),K2),kron(eye(r(3)),K3),kron(eye(r(4)),K4),kron(eye(r(5)),K5),kron(eye(r(6)),K6));
disp(KD)
disp(TC)
disp(TR)