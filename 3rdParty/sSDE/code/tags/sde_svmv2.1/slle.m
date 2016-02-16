% supervised LLE ALGORITHM based on "Supervised locally linear embedding"
%
% function [Y,WW,M] = sllev1(X,groupnames, K,alpha, tol,varargin)
% INPUT
% X = data as D x N matrix (D = dimensionality, N = #points)
% groupnames = the corresponding labels of the samples
% K = number of neighbors
% alpha = the supervision weight between two objective terms.
% d = embedding dimensionality
% tol = regularizer (defaults to 1e-4)
% OUTPUT
% Y = embedding as d x N matrix
% WW = the local weight matrix
% M = (I-WW)'(I-WW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y,WW,M] = slle(X,groupnames, K,d,alpha, tol,varargin)

% PAIRWISE DISTANCES
[D,N] = size(X);
X2 = sum(X.^2);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
% slle regularization
max_dist = max(distance(:));
[groupIndex, groupString] = grp2idx(groupnames);
Lam_mat = ones (length(groupIndex));
for i =1: length(groupString)
    Lam_mat(groupIndex==i, groupIndex==i) = 0;
end
Lam_mat = Lam_mat*alpha*max_dist +distance;
%; slle regularization 
% NEIGHBORS
for i=1:N
    distance(i,i)=max(distance(i,:));
end
[sorted,index] = sort(distance);
neighbors = index(1:K,:);

% RECONSTRUCTION WEIGHTS
if (nargin<6), tol=1e-4; end;
W = zeros(K,N);
WW=zeros(N);
for i=1:N
%   z = X(:,neighbors(:,i))-repmat(X(:,i),1,K);
%   C = z'*z;
  % slle
  z = Lam_mat(i, neighbors(:,i));
  C = 0.5* ( repmat(z', 1, K) + repmat(z, K, 1) - Lam_mat(neighbors(:,i), neighbors(:,i)) );
%   norm(C1 - C, 'fro')
  %; slle
  C = C + tol*trace(C)*eye(K)/K; % REGULARIZATION
  invC = inv(C);
  W(:,i) = sum(invC)'/sum(sum(invC));
  WW(neighbors(:,i),i)=W(:,i);
end;

% COST MATRIX
M = eye(N);
for i=1:N 
  w = W(:,i);
  j = neighbors(:,i);
  M(i,j) = M(i,j) - w';
  M(j,i) = M(j,i) - w;
  M(j,j) = M(j,j) + w*w';
end;

% CALCULATION OF EMBEDDING
options.disp = 0;
[Y,eigenvals] = eigs(M+1e-10*eye(N),d+1,0,options);
eigenvals = diag(eigenvals);
[eigenvals,indx] = sort(eigenvals);
eigenvals = eigenvals(2:d+1);
Y = Y(:,indx(2:d+1))'*sqrt(N);
Y=Y';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
