% Modified LLE ALGORITHM, ONLY computing the local linear weight
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% function [WW] = lle(X,K,dis)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% K = number of neighbors
% WW the local linear weight, each column contains the weight to reconstruct
% the corresponding sample. 
% dis the maximum distance that defines the neighborhood. NOT USED 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WW = llew(X,K,dis)

% PAIRWISE DISTANCES
[D,N] = size(X);
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

% NEIGHBORS
for i=1:N
    distance(i,i)=max(distance(i,:));
end
[sorted,index] = sort(distance);
neighbors = index(1:K,:);

% RECONSTRUCTION WEIGHTS
tol=1e-4; 
W = zeros(K,N);
WW=zeros(N,N);
for i=1:N
  z = X(:,neighbors(:,i))-repmat(X(:,i),1,K);
  C = z'*z;
  C = C + tol*trace(C)*eye(K)/K; % REGULARIZATION
  invC = inv(C);
  W(:,i) = sum(invC)'/sum(sum(invC));
  WW(neighbors(:,i),i)=W(:,i);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
