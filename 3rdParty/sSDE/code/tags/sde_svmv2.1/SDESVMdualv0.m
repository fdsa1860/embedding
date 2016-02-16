% SDE + soft margin SVM with y'*u==0 constraint
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUTs:
% St        the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% Sg        the high dimension test data, which is Ng-by-D, where Ng is the number
%           of samples
% NN        the number of local neightbors.
% Nl        the number of landmarks
% eps       the ratio of local geometric constraint relaxation [0 1].
% y         the label of each sample {-1 +1}.
% lam1      the weight between SDE and Maximum Margin objectives [0 1].
% lam2      the weight of the soft margin error in the Maximum Margin
%           objective.
% mode      the embedding mode. 'SS' indicates the semisupervised mode,
%           which is not fully tested.
% OUTPUT:
% K         the learned kernel/grammian matrix. 
% eat       the adjacent matrix used to solve the problem.
% indl      the index of the landmarks
% cvx_optval    the optimum computed by cvx
% cvx_status    the solver status for solving the SDP optimization
% t         the upper bound of the Maximum Margin objective.
% v1        the dual variables 
% v2        the dual variables 
% v3        the dual variable 
% cvx_result    each row records the results in each iteration, which is 
%               consisted of the number of active local geometric
%               constraints, the number of violated constraints ,the
%               objective function value, the summation of the amount of 
%               the violated constraints, the solving time until this 
%               iteration, the number of iterations used by the SDP solver.

function [K eta cvx_optval cvx_status t v1 v2 v3] = SDESVMdualv0(St,Sg, NN,eps,y,lam1, lam2, mode)
% clc
% initialization
time1= cputime;
K=[];  eta=[]; cvx_optval=[];  cvx_status=[]; 
t=[]; v1=[];  v2=[];  v3=[];
if lam1>1
    return;
end
if strcmp(mode,'SS')
    S1=double([St; Sg]);
else
    S1=double(St);
end
optvalcvx=0;
[N ~]=size(S1);
Nt=size(St,1);
Ng=size(Sg,1);
G=S1*S1';
relaxratio=eps;
%% Rank Minimization with Kernel Matrix
% Neighborhood Matrix
W=NNeighbour(S1',NN);
eta=W';
eta=double(abs(eta)>0);
[eta]=FullConn(S1, eta);
eta=(eta>0)|((eta'*eta)>0);
eta=eta|eta';
%%%%%%%%%%
[MC ,~]=FormConstraintAll(eta);
MGM=sum((MC*G).*MC,2);
D=diag(y);
ix_pos=find(y==1);
ix_neg=find(y==-1);
% % display('begin CVX')
%%
cvx_clear
cvx_begin sdp quiet
cvx_solver  sedumi %sdpt3 %sedumi  
cvx_precision medium
variable K(N,N) symmetric;
variable v1(Nt,1);
variable v2(Nt,1);
variable v3;
variable t;
minimize (((1-lam1)*t-lam1*trace(K))/sqrt((1-lam1)*lam1));%+n'*n
% minimize (((1-lam1)*t-lam1*(trace(K)-sum(sum(K(ix_pos,ix_pos)))/length(ix_pos) - sum(sum(K(ix_neg,ix_neg)))/length(ix_neg)))/sqrt((1-lam1)*lam1)) %+n'*n
subject to

sum(sum(K))==0;
% sum((MC*K).*MC,2)<MGM;
sum((MC*K).*MC,2)<= (1+eps)*MGM;
if eps <1
    sum((MC*K).*MC,2)>= (1-eps)*MGM;
end
% Constraint from SVM
if size(K,1)>Nt
    [K(1:Nt,1:Nt)        K(1:Nt,Nt+1:end)     D*(1+v1-v2+v3*y);...
        K(Nt+1:end,1:Nt)    K(Nt+1:end,Nt+1:end) zeros(N-Nt,1);...
        (1+v1-v2+v3*y)'*D'  zeros(1,N-Nt)        2*t-2*lam2*v2'*ones(Nt,1)]>=0;
    % K>=0;
else
    % [D*K(1:Nt,1:Nt)*D sqrt(2)/2*(1+v1-v2+v3*y); sqrt(2)/2*(1+v1-v2+v3*y)' t-lam2*v2'*ones(Nt,1)]>=0;
    [2*K(1:Nt,1:Nt) D*(1+v1-v2+v3*y);...
    (1+v1-v2+v3*y)'*D' t-lam2*v2'*ones(Nt,1)]==semidefinite(Nt+1);
    % [2*D*K(1:Nt,1:Nt)*D (1+v1-v2); (1+v1-v2)' t-lam2*v2'*ones(Nt,1)]>=0;
end
v1>=0;
v2>=0;
cvx_end
time2 = cputime;
display([datestr(clock)  ' runtime=  ' num2str(time2 -time1)]);
return;