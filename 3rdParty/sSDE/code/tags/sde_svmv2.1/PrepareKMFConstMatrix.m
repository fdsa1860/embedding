% function that compute the Kernel Factorization Matrix Q in the paper
% TODO: this function should be simplified. too many inputs.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       S1: data matrix, each column is a sample vector.
%       G: S1*S' 
%       W: the local linear weight matrix computed by "llew" function
%       eta: the local neiborhood connection matrix.
%       Nl: the number of landmarks
%       N: the number of all samples.
%       NN: the number of neighborhood.
% OUTPUT
%       Q: the factorized projection matrix.
%       QQ: Q'*Q
%       Qv: the matrix of (1'*Q)*(1'*Q)'.
%       S_etaQ: the local geometric constraints matrix times the
%       factorization matrix, S_eta*Q.
%       S_etaQv:  S_eta*Qv.
%       act_const_idx: the initialized  index of the active constraints.
%       dis: the vector of the input space local distance.
%       indl: the index of the landmarks. 
function [Q,QQ, Qv, S_etaQ, S_etaQv,act_const_idx, dis, indl ] = ...
    PrepareKMFConstMatrix(S1, G, W, eta,Nl, N, NN, varargin )
if length(varargin)>=1
    indl = varargin{1};
    [Q, indl] = ComputeFactorMatrix(W, Nl, N, NN, indl);
else
    min_val =Inf;
    Q_opt =[];
    for i =1:100
        [Q, indl] = ComputeFactorMatrix(W, Nl, N, NN);
        %     [Q, indl] = ComputeFactorMatrix(S1', eta, Nl, N, NN);
        Q(abs(Q)<1e-4)=0;
        Q=sparse(Q);
        if norm(G-Q*G(indl,indl)*Q','fro')/norm(G,'fro') < min_val
            min_val=norm(G-Q*G(indl,indl)*Q','fro')/norm(G,'fro');
            Q_opt = Q;
        end
        if min_val <0.001
            break;
        end
    end
    Q = Q_opt;
end

% [MC ,~]=FormConstraintAll(eta);
% MGM=sum((MC*G).*MC,2);
[S_eta, act_const_idx, act_const_idx2]=FormConstraint(eta, S1, indl);
dis=sum(S_eta*G.*S_eta,2);
% act_const_idx=act_const_idx| sum(abs(S_eta(:,indl)),2)>0;
S_etaQ=S_eta*Q;
QQ=Q'*Q;
Qv=sum(Q);
Qv=Qv'*Qv;
num_acc_const=sum(act_const_idx);

temp1= repmat(S_etaQ, [1, size(S_etaQ,2)]);
S_etaQv= repmat(S_etaQ, [size(S_etaQ,2),1]);
S_etaQv= reshape(S_etaQv, size(temp1));
S_etaQv=temp1.*S_etaQv;
% form a vectorize structure matrix
% idx=find(triu(ones(Nl)));
% temp=zeros(Nl);
% temp(idx)=1:length(idx);
% temp=temp';
% temp(idx)=1:length(idx);
% S_v=sparse(1:Nl^2,temp(:),ones(Nl^2,1),Nl^2,length(idx));
%
return;

% Compute Q according  "Nonlinear Dimensionality Reduction by Semidefinite 
% Programming and Kernel Matrix Factorization"
function [Q, indl] = ComputeFactorMatrix(W, Nl, N, NN, varargin)
if length(varargin)>=1
    indl = varargin{1};
else
    indl = randperm(N);
    indl = indl(1:Nl);
end
W=sparse(W');

indu=ones(1,N);
indu(indl)=0;
indu=find(indu==1);
phi=sparse(eye(N)-W);
phi=phi'*phi;
phi_uu=phi(indu,:);
phi_ul=phi_uu(:,indl);
phi_uu=phi_uu(:,indu);
Q=zeros(N,Nl);
Q(indl,:)=eye(Nl);
temp=inv(phi_uu+1e-15*eye(size(phi_uu)))*phi_ul;  % The equation in the paper lack a minus 
Q(indu,:)=-temp;
return

% function [Q, indl] = ComputeFactorMatrix(X, eta, Nl, N, NN)
% 
% dis = repmat(sum(X.^2, 1), [ size(X, 2), 1]);
% dis = dis + dis' - 2*X'*X;
% dis = sparse(eta.*dis);
% % select landmarks
% indl = randperm(N,Nl);
% indu= [1:N];
% indu(indl)=[];
% 
% for i=1:length(indl)
%     for j= 1:length(indu)
%         [ DIST(i,j), PATH(i,j) ] = graphkshortestpaths( dis, indu(j), indl(i),  1 );
%     end
% end
%     
% [sorted,index] = sort(DIST);
% neighbors = indl(index(1:NN,:));
% 
% % RECONSTRUCTION WEIGHTS
% tol=1e-4; 
% W = zeros(NN,N);
% WW=zeros(N,N);
% for i=1:N-Nl
%   z = X(:,neighbors(:,i))-repmat(X(:,indu(i)),1,NN);
%   C = z'*z;
%   C = C + tol*(1e-4+trace(C))*eye(NN)/NN; % REGULARIZATION
%   invC = inv(C);
%   W(:,indu(i)) = sum(invC)'/sum(sum(invC));
%   WW(neighbors(:,i),indu(i))=W(:,indu(i));
% end;
% Q(indl,:)=eye(Nl);
% Q(indu,:)= WW(indl,indu)';
% return