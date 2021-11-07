function [ K,W,J_all ] = algo2_synK(A,B2,B1,C,D,K,k,epsilon) 
%% matlab codes for Algorithm 2 in the manuscript "Structured controller synthesis through block-diagonal factorization and parameter space optimization"
%% YALMIP toolbox and the solver Gurobi are required to run this algorithm
%% authors: yuanliming@nimte.ac.cn, chensilu@nimte.ac.cn
%% set the iteration index
l = 0; 
%% run Algo1 to obtain the factorization results
[r,TC,TR] = algo1_facK(K,k); % k={K1,K2,...,KN}
%% transforamtion to decentralized feedback system (26)
bar_A = TR*(A)*((TR'*TR)^-1*TR'); 
bar_B2 = TR*B2*TC;
bar_B1 = TR*B1;
bar_C = C*((TR'*TR)^-1*TR'); 
bar_D = D*TC;
%% project the decentralized feedback system on a parameter space (29)
N = length(k);% the number of different sub-controllers
KD = [];
for i = 1:N
KD=blkdiag(KD,kron(eye(r(i)),k{i})); % structure of KD matrix eq.(14)
end
[bar_m,bar_n] = size(KD); % size of decentralized controller
F1 = [bar_A -bar_B2;zeros(bar_m,bar_n) zeros(bar_m,bar_m)];
Q1 = [bar_B1*bar_B1' zeros(bar_n,bar_m);
    zeros(bar_m,bar_n) zeros(bar_m,bar_m)]; 
R = blkdiag(bar_C'*bar_C,bar_D'*bar_D); 
%% define the paramer matrix, constraints and object fucntion
W1 = [];
W2 = [];
W3 = [];
for i=1:N % incorporate the structural constraints of KD in (13)
[mi,ni]=size(k{i}); 
W1 = blkdiag(W1,kron(eye(r(i)),sdpvar(ni,ni)));% Theorem 3, eq.(35)
W2 = blkdiag(W2,kron(eye(r(i)),sdpvar(ni,mi)));% Theorem 3, eq.(35)
W3 = blkdiag(W3,kron(eye(r(i)),sdpvar(mi,mi)));% Remark 6 (iii)
end
W = [W1  W2;
     W2' W3]; % construct the parameter matrix eq.(32) 
THETA1 = F1*W+W*F1'+Q1; % the matrix function for robust stabilization eq.(33)
%% initialization
P = [];
for i = 1:(bar_m+bar_n)
P = [P, W(i,i) >= 0];
end % initial constraints for a finite solution (positive diagonal element)
z = trace(R*W); % upper bound of H2 norm square eq.(40)
assign(W,0) % initial solution
ep1=0;
ep2=max((eig(Q1)));
%% cutting plane algorithm for serching the optimal W for Problem 2
J_all = [];
 while (-ep1) >= epsilon || (ep2) >= epsilon
    if ep1 <= -epsilon
    disp('ep1 violated') 
    v1 = mineigv(value(W));
    P = [P, v1'*W*v1 >=0]; % include the cutting plane according to violation of ep1
    end
    if ep2 >= epsilon
    disp('ep2 violated')
    v2 = maxeigv(value(THETA1(1:bar_n,1:bar_n))); 
    v2 = [v2;zeros(bar_m,1)];
    P = [P, v2'*THETA1*v2 <= 0]; % include the cutting plane according to violation of ep2
    end
  options = sdpsettings('solver','gurobi');
  optimize(P,z,options);% solve LP
  J = value(z);
  J_all = [J_all;J];% the evolution of the cost
  l=l+1;
  [ep1,~] = min((eig(value(W))));
  [ep2,~] = max((eig(value(THETA1(1:bar_n,1:bar_n)))));
 end
W=value(W); % the solution to Probnlem 2
%% Return the solution to Problem 1
K=TC*value(W2)'*(value(W1))^-1*TR; 
end
function [vector] = mineigv(W)
% obtain unit-norm eigenvector corresponding to the mimimum eigenvalue of a matrix
[~,pep1]=min(eig(W));
[Eig,~]=eig(W);
vector=Eig(:,pep1);%normlized eigenvector
end
function [vector] = maxeigv(W)
% obtain unit-norm eigenvector corresponding to the maximum eigenvalue of a matrix
[~,pep1]=max(eig(W));
[Eig,~]=eig(W);
vector=Eig(:,pep1);%normlized eigenvector
end
