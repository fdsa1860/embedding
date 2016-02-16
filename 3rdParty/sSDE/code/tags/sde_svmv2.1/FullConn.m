% Checking whether the graph is fully connected.
% If not, adding edges to connect the graph.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       X: data matrix, each row is a sample vector
%       eta: the neighborhood matrix. eta_ij =1 if and only if there is an
%           edge between i-th and j-th sample.
% OUTPUT
%       eta: the output neighborhood matrix in which the whole graph is 
%           fully connected. eta_ij =1 if and only if there is an
%           edge between i-th and j-th sample.
function [eta]=FullConn(X, eta)
eta=eta-diag(diag(eta));
S=1:size(X,1);
k=1;
% finding all the clips via depth first traverse 
while (~isempty(S))
    label=zeros(size(S));
    label=color(1,label,eta(S,S)); % label the clips with a index
    idx{k}=S(find(label)); % extract the index of samples in the clip.
    group{k}=X(idx{k},:);
    
    S(find(label))=[];
    k=k+1;
end
if length(group)==1
    return;
end
% compute the distance between group
dis=inf(length(group));
idxpair=cell(length(group),length(group));
for i=1:length(group)
    for j=i+1:length(group)
        temp=repmat(sum(group{i}.^2,2),1,size(group{j},1));
        temp=temp+ repmat(sum(group{j}.^2,2)',size(group{i},1),1);
        temp=temp-2*group{i}*group{j}';
        dis(i,j)=min(temp(:));
        dis(j,i)=min(temp(:));
        [idxi,idxj]=find(temp==min(temp(:)));
        idxpair{i,j}(1)=idx{i}(idxi(1));
        idxpair{i,j}(2)=idx{j}(idxj(1));
        idxpair{j,i}(1)=idx{j}(idxj(1));
        idxpair{j,i}(2)=idx{i}(idxi(1));
    end
end


% merge the groups
for k=length(group):-1:2
    ix=1:k;
    [idxi, idxj]=find(dis==min(dis(:)));
    idxi=idxi(1);idxj=idxj(1);
    eta(idxpair{idxi,idxj}(1),idxpair{idxi,idxj}(2))=1;
    eta(idxpair{idxi,idxj}(2),idxpair{idxi,idxj}(1))=1;
    if k<=2;
        return;
    end
    ixs=[idxi  idxj];
    ix(ixs)=[];
    
    temp=dis(ix,:);
    [temp, ixtemp]=min(temp(:,ixs),[],2);
    dis=[dis(ix,ix) temp; temp' inf];
    newidxpair=cell(k-1,k-1);
    newidxpair(1:k-2,1:k-2)=idxpair(ix,ix);
    for i=1:k-2
        newidxpair{i,k-1}=idxpair{ix(i),ixs(ixtemp(i))};
        newidxpair{k-1,i}=idxpair{ixs(ixtemp(i)),ix(i)};
    end
    newidxpair{end,end}=[nan nan];
    idxpair=newidxpair;
end
return;
% Depth First Traverse
function [label]=color(i,label,eta)
if label(i)==1
    return;
end
label(i)=1;
idx=find(eta(i,:));
for k=1:length(idx)
    if label(idx(k))~=1;
        label=color(idx(k),label,eta);
    end
end


