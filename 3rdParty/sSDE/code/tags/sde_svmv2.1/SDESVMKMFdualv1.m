% sSDE_l function
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
function [K eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdualv1(St,Sg, NN,Nl,eps,y,lam1, lam2, mode,varargin)
% clc
% initialization
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
% compute the SVs
% idx_SV = FindSVs(St, y, lam2/4);
% Nl =length(idx_SV);

if length(varargin)>=1
    W=llew(S1',NN,varargin{1});
else
    W=llew(S1',NN);
end

eta=double(abs(W)>0)';
eta=double(abs(eta)>0);
% make sure of the the samples are all connected. Otherwise the problem is
% unbounded.
[eta]=FullConn(S1, eta); 
eta=(eta>0)|((eta'*eta)>0);
eta=eta|eta';   % The definition of eta matrix in SDE paper.

% initialize the subset of the local isometric inequality constraints. And
% compute the factorization matrix Q.
[Q,QQ, Qv, S_etaQ, S_etaQv,act_const_idx, dis , indl] = PrepareKMFConstMatrix(S1, G, W, eta,Nl, N, NN );

D=diag(y);
ix_pos=find(y==1);  % index of the positive samples
ix_neg=find(y==-1); % index of the negative samples
relax_zero = 1e0;   % relaxed value for sum(K(:))==0 constraint
%%
num_act_const_step=200;
cvx_result = [];
cpu_time=0;
cvx_clear;
while 1
    % iteratively solving the SDP with part of the local isometric
    % inequality constraints
    disp('CVX begin')
    cvx_begin sdp quiet
    cvx_solver  sedumi%sdpt3%
    cvx_precision medium
    variable L(Nl,Nl) symmetric;
    variable v1(Nt,1);
    variable v2(Nt,1);
    variable v3;
    variable t;
    variable v(Nl, 1);
    
    minimize( ( (1-lam1)*t-lam1*trace(L*QQ) )/sqrt((1-lam1)*lam1) );
    subject to
    sum(L(:).*Qv(:))<=relax_zero; % relaxing the constraint: sum(L(:).*Qv(:))==0;
    S_etaQv(act_const_idx,:)*L(:)<=(1+eps) * dis(act_const_idx);
    [2*L  v; v'  t-lam2*v2'*ones(N,1)]==semidefinite(Nl+1);
    v1>=0;
    v2>=0;
	D*Q(1:Nt,:)*v == 1 + v1 - v2 + v3*y;
    cvx_end
    disp(cvx_status)
    %%
    if (~(strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) || sum(isnan(L(:)))>0 )
        % if cvx fail to solve the problem, reinitialize the constriants
        % and solve it again.
        [Q,QQ, Qv, S_etaQ, S_etaQv,act_const_idx, dis , indl] = PrepareKMFConstMatrix(S1, G, W, eta,Nl, N, NN );
        cvx_result = [];
        continue;
    end

    diff=S_etaQv*L(:)-dis-1e-4;
    ix_g0 =find(diff>0 & act_const_idx ==0);
    if isempty(ix_g0)
        break;
    end
    % find the most violated local isometric constraints.
    diff = diff(ix_g0);
    [~,ix]=sort(diff,'descend');
    cpu_time = cpu_time + cvx_cputime;
    cvx_result = [cvx_result; sum(act_const_idx) length(ix) cvx_optval full(sum(diff)) cpu_time cvx_slvitr];
%     display(cvx_result(end,:))
    if length(ix)>num_act_const_step
        act_const_idx(ix_g0(ix(1:num_act_const_step)))=1;
    else
        act_const_idx(ix)=1;
    end
    if (size(cvx_result, 1)>1) && (cvx_result(end,3) - cvx_result(end-1,3))< 1e-4 * abs(cvx_result(end-1,3))
        break;
    end
     clear L v1 v2 v3 t v;
end
K =Q*L*Q';
display('FINISHED')
return;
