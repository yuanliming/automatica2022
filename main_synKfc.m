A=zeros(5,5);
B1=eye(5);
B2=eye(5);
C=[1,  0,  0,  0,  0; 
   1, -1,  0,  0,  0;
   0,  1, -1,  0,  0; 
   0,  0,  1, -1,  0; 
   0,  0,  0,  1, -1; 
   0,  0,  0,  0,  1;
    zeros(5,5)];
D=[zeros(6,5);eye(5)];
syms k11 k12 k23 k34 k45 k55 
Kfc=[k11+k12        -k12           0           0          0;
        -k12     k12+k23        -k23           0          0;
           0        -k23     k23+k34        -k34          0;
           0           0        -k34     k34+k45       -k45;
           0           0           0        -k45    k45+k55;];
K1=k11;K2=k12;K3=k23;K4=k34;K5=k45;K6=k55;
epsilon=1*10^-6;
[Kfc, W, J]=algo2_synK(A,B2,B1,C,D,Kfc,{K1,K2,K3,K4,K5,K6},epsilon);

disp(W)
disp(Kfc)

plot(J,'linewidth',2)
xlabel('Iteration number','fontsize',16);
ylabel('tr(RW)','fontsize',16,'Interpreter','latex');
set(gca,'FontSize',16) 
set(gca,'LineWid',0.8)
set(gca,'xcolor',[0 0 0])
set(gca,'ycolor',[0 0 0])
grid on