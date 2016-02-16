%% out-of-sample extension with RBF approach
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
% St        the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% LD_test   the training data representatives in the output space, each
%           column is a embedded data sample. 
% Sg        the high dimension test data, which is Ng-by-D, where Ng is the number
%           of samples
% y         the label of each sample {-1 +1}.
% nCenters  the number of clusters,parameter for learning RBF mapping.
% LLM       the local linear projection matrices, for the RBF mapping.
% kmeans_space  the centers of each cluster, for RBF mapping.
% OUTPUT
% LD_test   the testing data representatives in the output space, each
%           column is a embedded data sample.
% cluster_label     the label of cluster for each sample
% LLM       the structure that store the local cluster information.
function [LD_test, cluster_label,LLM] = ...
    NLMap(St,LD_train,Sg, y,nCenters, LLM,kmeans_space)
LD_test=[]; cluster_label=[];
%% RBF mapping
% nCenters=40;
% c=0.1;
% RBF = trainGRBF(St',LD_train',nCenters,c);
% Y=testGRBF(RBF,Sg',c);
% LD_test=Y';
%%
if isempty(LLM)
    if strcmp(kmeans_space,'HD')
        [LLM, idx_pos, idx_neg]=LLM_trainv1(St,LD_train,nCenters,y);
    elseif strcmp(kmeans_space,'LD')
        [LLM, idx_pos, idx_neg]=LLM_trainv2(St,LD_train,nCenters,y);
    end
end
[LD_test,cluster_label]=LLM_test(Sg,LLM);
%% lle mapping
% [WW]= lleweight(St, Sg,NN);
% LD_test=WW*LD_train;
return;
% DEPRECATED
% function [WW]= lleweight(S_train, S_test,NN)
% WW=zeros(size(S_train,1), size(S_test,1));
% dis = repmat(sum(S_train.^2,2),1,size(S_test,1));
% dis = dis+ repmat(sum(S_test.^2,2)',size(S_train,1),1);
% dis = dis - 2*S_train*S_test';
% 
% [sorted, index] = sort(dis,1);
% neighbors = index(1:NN,:);
% 
% for i=1:size(S_test,1)
%     z =S_train(neighbors(:,i),:)-repmat(S_train(i,:),NN,1);
%     C = z*z';
%     C = C + 1e-10*trace(C)*eye(NN)/NN; % REGULARIZATION
%     invC = inv(C);
%     W(:,i) = sum(invC)'/sum(sum(invC));
%     WW(neighbors(:,i),i)=W(:,i);
% end
% WW =WW';
% return;

function [LLM, idx_pos, idx_neg]=LLM_trainv1(St,LD_train,nCenters,y)
%%
St_pos=St;
St_neg=[];
idx_neg =[];
centersX_neg=[];
nCenters_pos=nCenters;
nCenters_neg=0;
if ~isempty(y)
    St_pos(y==-1,:)=[];
    St_neg=St(y==-1,:);
    nCenters_pos=max(round(nCenters*sum(y==1)/length(y)),1);
    nCenters_neg=nCenters-nCenters_pos;
    if nCenters_neg==0
        nCenters_neg=1;
        nCenters_pos=nCenters_pos-1;
    end
end
[ idx_pos centersX_pos sumd] = rep_litekmeans( St_pos, nCenters_pos,50);

if ~isempty(St_neg)
    [ idx_neg centersX_neg sumd] = rep_litekmeans( St_neg, nCenters_neg,50);
end
idx=zeros(size(y));
idx(y==1)=idx_pos;
idx(y==-1)=idx_neg + nCenters_pos;
centersX=[centersX_pos; centersX_neg];

%%
P=cell(nCenters,1);
for i=1:nCenters
    ix=find(idx==i);
    S=St(ix,:)-repmat(centersX(i,:),length(ix),1);
    s=LD_train(ix,:);
    centersY(i,:)=mean(s,1);
    s=s-repmat(centersY(i,:),length(ix),1);
    P{i}=(S+1e-11*eye(size(S)))\(s);
end
LLM.nCenters=nCenters;
LLM.centersX=centersX;
LLM.centersY=centersY;
LLM.P=P;
%%
return;
function [LLM, idx_pos, idx_neg]=LLM_trainv2(St,LD_train,nCenters,y)
%%
St_pos=LD_train;
St_neg=[];
idx_neg =[];
centersX_neg=[];
nCenters_pos=nCenters;
nCenters_neg=0;
if ~isempty(y)
    St_pos(y==-1,:)=[];
    St_neg=LD_train(y==-1,:);
    nCenters_pos=max(round(nCenters*sum(y==1)/length(y)),1);
    nCenters_neg=nCenters-nCenters_pos;
    if nCenters_neg==0
        nCenters_neg=1;
        nCenters_pos=nCenters_pos-1;
    end
end

[ idx_pos centersX_pos sumd] = rep_litekmeans( St_pos, nCenters_pos, 50);

if ~isempty(St_neg)
    [ idx_neg centersX_neg sumd] = rep_litekmeans( St_neg, nCenters_neg, 50);
end
idx=zeros(size(y));
idx(y==1)=idx_pos;
idx(y==-1)=int16(idx_neg) + int16(nCenters_pos);
for i=1:max(idx)
    centersX(i,:)=mean(St(idx==i,:),1);
end

%%
nCenters=size(centersX,1);
P=cell(nCenters,1);
for i=1:nCenters
    ix=find(idx==i);
    S=St(ix,:)-repmat(centersX(i,:),length(ix),1);
    s=LD_train(ix,:);
    centersY(i,:)=mean(s,1);
    s=s-repmat(centersY(i,:),length(ix),1);
    P{i}=(S+1e-11*eye(size(S)))\(s);
end
LLM.nCenters=nCenters;
LLM.centersX=centersX;
LLM.centersY=centersY;
LLM.P=P;
%%
return;
function [s,cluster_label]=LLM_test(X,LLM)
dis= repmat(sum(X.^2,2),1, LLM.nCenters);
dis= dis + repmat(sum(LLM.centersX.^2,2)',size(X,1), 1);
dis= dis - 2*X*LLM.centersX';
[temp ,ix]=sort(dis');
for i=1:size(X,1)
    s(i,:) = (X(i,:) - LLM.centersX(ix(1,i),:))*LLM.P{ix(1,i)} + LLM.centersY(ix(1,i),:);
end
cluster_label=ix(1,:);
return;

function [idx, centersX, sumd]= rep_litekmeans(X, k, it)
X=X';
sumd=Inf;
centersX=[];
idx=[];
n = size(X,2);

for i=1:it
    label= litekmeans(X, k);
    [u,~,label] = unique(label);   % remove empty clusters
    kk = length(u);
    E = sparse(1:n,label,1,n,kk,n);  % transform label into indicator matrix
    m = X*(E*spdiags(1./sum(E,1)',0,kk,kk));    % compute m of each cluster
    [dis,label] = min(bsxfun(@minus,dot(m,m,1)'/2,m'*X),[],1); % assign samples to the nearest centers
    dis= sum(dis);
    if dis<sumd
        idx=label';
        centersX=m';
        sumd=dis;
    end
end

return;