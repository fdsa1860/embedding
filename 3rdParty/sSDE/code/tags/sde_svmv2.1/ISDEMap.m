%% out-of-sample extension with local linear assumption
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
% NN        the number of local neightbors.
% OUTPUT
% LD_test   the testing data representatives in the output space, each
%           column is a embedded data sample.
function [LD_test] = ISDEMap(St,LD_train, Sg, NN, varargin)
LD_test=[]; 
if ~isempty(varargin) && strcmp(varargin{1},'LLW')
    [WW]= lleweight(St, Sg,NN);
    LD_test=WW*LD_train;
    return;
end
% initialization

S1=double([St; Sg]);
[N ~]=size(S1);
X=S1-kron(ones(N,1),mean(S1,1));
G=X*X';
X=X./sqrt(trace(G)/N);
G=G/trace(G)*N;
optvalcvx=0;

Nt=size(St,1);
Ng=size(Sg,1);
%% Rank Minimization with Kernel Matrix
W=llew(S1',NN);

indl=1:Nt;
W=sparse(W');
indu=1+Nt:Nt+Ng;
phi=sparse(eye(N)-W);
phi=phi'*phi;
phi_uu=phi(indu,:);
phi_ul=phi_uu(:,indl);
phi_uu=phi_uu(:,indu);
Q=zeros(N,Nt);
Q(indl,:)=eye(Nt);
temp=inv(phi_uu+1e-15*eye(size(phi_uu)))*phi_ul;  % The equation in the paper lack a minus 
% Q(indu,:)=-temp;
LD_test=(-temp*LD_train);
return;

function [WW]= lleweight(S_train, S_test,NN)
WW=zeros(size(S_train,1), size(S_test,1));
dis = repmat(sum(S_train.^2,2),1,size(S_test,1));
dis = dis+ repmat(sum(S_test.^2,2)',size(S_train,1),1);
dis = dis - 2*S_train*S_test';
dis(dis ==0)=Inf;
[sorted, index] = sort(dis,1);
neighbors = index(1:NN,:);
if (nargin<4), tol=1e-4; end;
for i=1:size(S_test,1)
    z =S_train(neighbors(:,i),:)-repmat(S_test(i,:),NN,1);
    C = z*z';
    C = C + tol*trace(C)*eye(NN)/NN; % REGULARIZATION
    invC = inv(C);
    W(:,i) = sum(invC)'/sum(sum(invC));
    WW(neighbors(:,i),i)=W(:,i);
end
WW =WW';
return;