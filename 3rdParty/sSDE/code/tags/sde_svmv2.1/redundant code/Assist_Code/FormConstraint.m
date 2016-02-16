% function that form the local geometric constraint.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       S1: data matrix, each column is a sample vector.
%       eta: the local neiborhood connection matrix.
%       indl: the index of the landmark samples.
% OUTPUT
%       S_eta: coefficient vectors for the local isometric constraints. So
%           that for each edge (eta_ij =1) there is a row in S_eta that the
%           i-th and j-th entry in that row are 1 and -1 respectively and
%           all the other entries are 0.
%       acc_idx1: the index of the active local geometric constraints.

function [S_eta, acc_idx1,acc_idx2]=FormConstraint(eta, S1, indl)
    acc_idx2=[];
    eta=sparse(eta);
    eta=eta|eta';
    eta=eta-diag(diag(eta));
    eta=triu(eta);
%     for i=1:length(indl)
%         for j= i:length(indl)
%             [ DIST(i,j), PATH(i,j) ] = graphkshortestpaths( eta, indl(i), indl(j), 1 );
%             DIST(j,i)= DIST(i,j);
%             PATH(j,i)=PATH(i,j);
%         end
%     end
%     [~, ix]=min(sum(DIST));
%     PATH=PATH(ix,:);
%     acc_eta=zeros(size(eta));
%     for i=1:length(indl)
%         for j=1:length(PATH{i})-1
%             acc_eta(PATH{i}(j),PATH{i}(j+1))=1;
%         end
%     end
    S= sum(S1.^2,2);
    S =repmat(S, [1 size(S1,1)]);
    S = S+ S'-2*S1*S1';
    S =eta.*S;
    BGobj = biograph(S);
    [Tree, pred] = minspantree(BGobj);
    acc_eta =Tree>0;
    [acc_eta] = PruneMinTree(acc_eta, indl);
    
    acc_eta=acc_eta|acc_eta';
    acc_idx1=acc_eta(find(eta(:)));
    
    N_eta=sum(eta(:));
    [idxx idxy]=ind2sub(size(eta),find(eta(:)));
    temp=sparse([1:length(idxx)],idxx,ones(size(idxx)),N_eta, size(eta,1));
    S_eta=temp+sparse([1:length(idxy)],idxy,-ones(size(idxy)),N_eta, size(eta,1));
    
%     for i=1:length(indl)
%         [dist(i,:), path(i,:), ~]=shortestpath(bg2,indl(i),indl);
%     end
%     [~, ix]=min(sum(dist));
%     path=path(ix,:);
%     acc_eta=zeros(size(eta));
%     for i=1:length(indl)
%         acc_eta(indl(ix),path{i}(1))=1;
%         for j=1:length(path{i})-1
%             acc_eta(path{i}(j),path{i}(j+1))=1;
%         end
%     end
%     acc_eta=acc_eta|acc_eta';
%     acc_idx1=acc_eta(find(eta(:)));
%     [Tree, pred] = minspantree(bg2);
%     T=Tree+Tree';
%     acc_idx2=T(find(eta(:)));
return;
% eliminate the edge between non landmark node and keep the tree connected
function [eta] = PruneMinTree(eta, indl)
idxl = zeros(size(eta,1),1);
idxl(indl) = 1;
while 1
    idx_leaf = (sum(eta, 2) == 1) & (idxl==0);
    if sum(idx_leaf) ==0
        break;
    end
    eta(idx_leaf,:) = 0;
end
return;