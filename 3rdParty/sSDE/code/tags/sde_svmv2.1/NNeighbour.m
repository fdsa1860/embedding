% Find the K nearest neighbors for each sample
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       X: data matrix, each row is a sample vector
%       K: the size of neighborhood.
%       distance (optional): the predefined maximum distance of the neighborhood. 
% OUTPUT
%       WW: the output neighborhood matrix in which the whole graph is 
%           fully connected. WW_ij =1 if and only if the j-th sample is the
%           neighbor of the i-th sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WW = NNeighbour(X,K, varargin)

% PAIRWISE DISTANCES
[D,N] = size(X);
WW=zeros(N,N);
optargin = size(varargin,2);
if optargin>1
    warning('too many input arguments');
    return;
end
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

% NEIGHBORS
for i=1:N
    distance(i,i)=max(distance(i,:));
end
[sorted,index] = sort(distance);
neighborhood=cell(1,N);
tol=1e-4;
for i=1:N
    KK=K;
    if optargin==1 && sorted(KK,i)>varargin{1}
        KK=sum(double(sorted(:,i)<=varargin{1}));
    end
    neighbors{i} = index(1:KK,i);
    WW(neighbors{i},i)= 1;
end;
% WW(2:end,1:end-1)=WW(2:end,1:end-1)+eye(N-1);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
