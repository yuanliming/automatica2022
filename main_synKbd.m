A= [0.1975    0.7375    0.8957    0.5437    0.7883    0.4944;
    0.8620    0.3033    0.2741    0.3867    0.8529    0.0312;
    0.1256    0.0430    0.9979    0.8222    0.4856    0.8276;
    0.6456    0.8355    0.8344    0.5953    0.8738    0.2643;
    0.4373    0.3693    0.7913    0.7816    0.3338    0.6775;
    0.6126    0.6613    0.6566    0.9877    0.2053    0.7939];
B1=eye(6);
B2=[0.6818    0.2091    0.0875    0.5936;
    0.6508    0.2725    0.3089    0.0438;
    0.2373    0.7758    0.2309    0.4249;
    0.4774    0.3314    0.9092    0.5216;
    0.9364    0.6031    0.9369    0.8403;
    0.2411    0.1840    0.0319    0.6250];
C=[eye(6);
    zeros(4,6)];
D=[zeros(6,4);
    eye(4)];
syms k1 k2 k3 k4 k5 k6 k7 k8 k9 
Kbd=[k1 k2 k3 k1 k2 k3;
     k4 k5 k6 k4 k5 k6;
     0  0  0  k1 k2 k3;
     0  0  0  k4 k5 k6;];
K1=[k1 k2 k3;k4 k5 k6];
epsilon=1*10^-6;
[Kbd, W, J]=algo2_synK(A,B2,B1,C,D,Kbd,{K1},epsilon);

disp(W)
disp(Kbd)

plot(J,'linewidth',2)
xlabel('Iteration number','fontsize',16);
ylabel('tr(RW)','fontsize',16,'Interpreter','latex');
set(gca,'FontSize',16) 
set(gca,'LineWid',0.8)
set(gca,'xcolor',[0 0 0])
set(gca,'ycolor',[0 0 0])
grid on