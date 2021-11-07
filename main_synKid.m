m1=15.28; % the mass of carriage 1
m2=15.28; % the mass of carriage 2
mL=22.2; % the mass of the crossbeam
JL=5.8; % the moment of inertia of the crossbeam
L=1.514; % the length of the crossbeam
f1=8; % coulomb frictions in carriage 1
f2=12; % coulomb frictions in carriage 2
KF=85.5; % for constant of the linear motors
M=[m1+mL/4+JL/L^2    mL/4-JL/L^2;
   mL/4-JL/L^2       m2+mL/4+JL/L^2];
A=[0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   zeros(2,6)];
B1=[0 ; 
    0;
    0 ;
    0 ;
    M^-1*[-f1;
    -f2]];
B2=[zeros(2,2) zeros(2,1);
    zeros(2,2) zeros(2,1);
    M^-1*[KF 0    -1/L;
          0    KF  1/L];
    ];
C=[eye(2) zeros(2,4);
    zeros(3,6)];
D=[zeros(2,3);diag([0.1,0.1,0.001])];
syms ki1 kp1 kd1 ki2 kp2 kd2 twokm_L
Kid=[ki1 0  kp1  0  kd1 0;
   0  ki2 0   kp2 0  kd2;
   0  0  -twokm_L twokm_L 0  0;]; 
K1=[ki1 kp1 kd1];       
K2=[ki2 kp2 kd2];       
K3=twokm_L;               
epsilon=1*10^-6;
[Kid, W, J]=algo2_synK(A,B2,B1,C,D,Kid,{K1,K2,K3},epsilon);

disp(W)
disp(Kid)

plot(J,'linewidth',2)
xlabel('Iteration number','fontsize',16);
ylabel('tr(RW)','fontsize',16,'Interpreter','latex');
set(gca,'FontSize',16)
set(gca,'LineWid',0.8)
set(gca,'xcolor',[0 0 0])
set(gca,'ycolor',[0 0 0])
grid on
