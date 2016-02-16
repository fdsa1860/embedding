% SDE + soft margin SVM with y'*u==0 constraint
% INPUTs:
% St        the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% Sg        the high dimension test data, which is Ng-by-D, where Ng is the number
%           of samples
% NN        the number of local neightbors.

function [K eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdualv0(St,Sg, NN,Nl,eps,y,lam1, lam2, mode,varargin)
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

% compute factorization 
if length(varargin)>=1
    W=llew(S1',NN,varargin{1});
else
    W=llew(S1',NN);
end

eta=double(abs(W)>0)';
eta=double(abs(eta)>0);
[eta]=FullConn(S1, eta);
eta=(eta>0)|((eta'*eta)>0);
eta=eta|eta';


if Nl>N*0.9
    [S_eta ,~]=FormConstraintAll(eta);
    Q = sparse(eye(N));
    S_etaQ=S_eta;
    QQ=Q'*Q;
    Qv=sum(Q);
    Qv=Qv'*Qv;
    
    temp1= repmat(S_etaQ, [1, size(S_etaQ,2)]);
    S_etaQv= repmat(S_etaQ, [size(S_etaQ,2),1]);
    S_etaQv= reshape(S_etaQv, size(temp1));
    S_etaQv=temp1.*S_etaQv;
    
    act_const_idx =ones(size(S_eta,1),1);
    dis = S_etaQv*G(:);
    indl = 1:N;
else
    
    if ~isempty(varargin)
        [Q,QQ, Qv, S_etaQ, S_etaQv,act_const_idx, dis , indl] = PrepareKMFConstMatrix(S1, G, W, eta,Nl, N, NN, varargin{1});
    else
        [Q,QQ, Qv, S_etaQ, S_etaQv,act_const_idx, dis , indl] = PrepareKMFConstMatrix(S1, G, W, eta,Nl, N, NN);
    end
end



D=diag(y);
ix_pos=find(y==1);
ix_neg=find(y==-1);
relax_zero = 1e0;

cvx_result = [];
% compute the SVs
% idx_SV = FindSVs(St, y, lam2/4);
% Nl =length(idx_SV);cpu_time=0;
%%
cvx_clear;
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
% SDE constraint
sum(L(:).*Qv(:)) <= relax_zero;
S_etaQv*L(:)<=(1+eps) * dis;
% svm constraint
[2*L  v; v'  t-lam2*v2'*ones(N,1)]==semidefinite(Nl+1);
v1>=0;
v2>=0;
D*Q(1:Nt,:)*v == 1 + v1 - v2 + v3*y;
cvx_end

cvx_cputime
cpu_time = cvx_cputime;
cvx_result = [cvx_result; cvx_optval cpu_time ];

K =Q*L*Q';
return;
% 
% %%
% function [idx] = FindSVs(Training, y, C)
% 
% D=diag(y);
% class=y; class(class==-1)=2;
% u=zeros(size(y));
% % train SVM
% options = statset('MaxIter',1e7);
% SVMStruct = svmtrain(Training,class,'kernel_function','polynomial',...
%     'boxconstraint',C, 'autoscale', false,'options',options);
% u(SVMStruct.SupportVectorIndices,1)=SVMStruct.Alpha./y(SVMStruct.SupportVectorIndices);
% gamma=-SVMStruct.Bias;
% e=1-D*(Training*Training'*D*u-gamma);
% 
% idx = SVMStruct.SupportVectorIndices;
% return;
% 
% function [L, v1, v2, v3, t] = solver_cvx0( Nl, Nsv, y, dis, inv_DQ_sv, idx_SV,  Qv,  QQ, S_etaQ, act_const_idx, lam1, lam2)
% %%
% disp('CVX begin')
% %     SDE MF algorithm
% %     cvx_begin sdp
% %     cvx_solver  sedumi%sdpt3%
% %     cvx_precision medium
% %     variable L(Nl,Nl) symmetric;
% %
% %     minimize(-trace(L*QQ));
% %     subject to
% %     sum(L(:).*Qv(:))<1e1;
% %     L==semidefinite(Nl);
% %     sum(S_etaQ(act_const_idx,:)*L.*S_etaQ(act_const_idx,:),2)<=dis(act_const_idx);
% %     cvx_end
% %     disp('cvx end')
% 
% cvx_begin sdp
% cvx_solver  sedumi%sdpt3%
% cvx_precision medium
% variable L(Nl,Nl) symmetric;
% variable v1(Nsv,1);
% variable v2(Nsv,1);
% variable v3;
% variable t;
% %     expression L_SV(Nsv, Nsv);
% %     L_SV = DQ_sv*L*DQ_sv';
% minimize( ( (1-lam1)*t-lam1*trace(L*QQ) )*1e3 );
% subject to
% sum(L(:).*Qv(:))<relax_zero;
% sum(S_etaQ(act_const_idx,:)*L.*S_etaQ(act_const_idx,:),2)<dis(act_const_idx);
% % [2*L Dl*(1+v1-v2+v3*y); (1+v1-v2+v3*y)'*Dl t-lam2*v2'*ones(Nt,1)]>0;
% %     [2*L_SV (1+v1-v2+v3*y(idx_SV));
% %         (1+v1-v2+v3*y(idx_SV))' t-lam2*v2'*ones(Nsv,1)]==semidefinite(Nsv+1);
% [2*L inv_DQ_sv*(1+v1-v2+v3*y(idx_SV));
%     (1+v1-v2+v3*y(idx_SV))'*inv_DQ_sv' t-lam2*v2'*ones(Nsv,1)]==semidefinite(Nsv+1);
% v1>=0;
% v2>=0;
% cvx_end
% return;
% 
% function [L, v1, v2, v3, t] = solver_cvx1( Nl, N, y, dis,  Qv,  Q, S_etaQ, act_const_idx, lam1, lam2)
% %%
% disp('CVX begin')
% cvx_begin sdp
% cvx_solver  sedumi%sdpt3%
% cvx_precision medium
% variable L(Nl,Nl) symmetric;
% variable v1(N,1);
% variable v2(N,1);
% variable v3;
% variable t;
% variable alpha(Nl, 1);
% 
% minimize( ( (1-lam1)*t-lam1*trace(L*QQ) )*1e3 );
% subject to
% sum(L(:).*Qv(:))<relax_zero;
% sum(S_etaQ(act_const_idx,:)*L.*S_etaQ(act_const_idx,:),2)<dis(act_const_idx);
% % [2*L Dl*(1+v1-v2+v3*y); (1+v1-v2+v3*y)'*Dl t-lam2*v2'*ones(Nt,1)]>0;
% %     [2*L_SV (1+v1-v2+v3*y(idx_SV));
% %         (1+v1-v2+v3*y(idx_SV))' t-lam2*v2'*ones(Nsv,1)]==semidefinite(Nsv+1);
% [2*L alpha;
%     alpha' t-lam2*v2'*ones(N,1)]==semidefinite(N+1);
% v1>=0;
% v2>=0;
% Q*alpha == 1 + v1 - v2 + v3*y;
% cvx_end
% return;
% 
% % function [L, v1, v2, v3, t] = solver_sedumi( Nl, Nsv, y, dis, inv_DQ_sv, idx_SV,  Qv,  QQ, S_etaQ, act_const_idx, lam1, lam2)
% % %%
% % disp('sedumi begin')
% % %     SDE MF algorithm
% % %     cvx_begin sdp
% % %     cvx_solver  sedumi%sdpt3%
% % %     cvx_precision medium
% % %     variable L(Nl,Nl) symmetric;
% % %
% % %     minimize(-trace(L*QQ));
% % %     subject to
% % %     sum(L(:).*Qv(:))<1e1;
% % %     L==semidefinite(Nl);
% % %     sum(S_etaQ(act_const_idx,:)*L.*S_etaQ(act_const_idx,:),2)<=dis(act_const_idx);
% % %     cvx_end
% % %     disp('cvx end')
% % 
% % prob=sdmpb('SDSVMKMF sedumi solver');
% % [prob, idx_L] = sdmvar(prob, Nl, 's', 'L');
% % [prob, idx_v1] = sdmvar(prob, Nsv, 1, 'v1');
% % [prob, idx_v2] = sdmvar(prob, Nsv, 1, 'v2');
% % [prob, idx_v3] = sdmvar(prob, 1, 1, 'v3');
% % [prob, idx_t] = sdmvar(prob, 1, 1, 't');
% % 
% % %     expression L_SV(Nsv, Nsv);
% % %     L_SV = DQ_sv*L*DQ_sv';
% % minimize( ( (1-lam1)*t-lam1*trace(L*QQ) )*1e3 );
% % subject to
% % sum(L(:).*Qv(:))<relax_zero;
% % sum(S_etaQ(act_const_idx,:)*L.*S_etaQ(act_const_idx,:),2)<dis(act_const_idx);
% % % [2*L Dl*(1+v1-v2+v3*y); (1+v1-v2+v3*y)'*Dl t-lam2*v2'*ones(Nt,1)]>0;
% % %     [2*L_SV (1+v1-v2+v3*y(idx_SV));
% % %         (1+v1-v2+v3*y(idx_SV))' t-lam2*v2'*ones(Nsv,1)]==semidefinite(Nsv+1);
% % [2*L inv_DQ_sv*(1+v1-v2+v3*y(idx_SV));
% %     (1+v1-v2+v3*y(idx_SV))'*inv_DQ_sv' t-lam2*v2'*ones(Nsv,1)]==semidefinite(Nsv+1);
% % v1>=0;
% % v2>=0;
% % return;