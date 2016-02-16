% %% calculate the dimensionality of the kernel matrix
% function [LD_train, dim] = CalDim(K_train, dim_thresh)
% [U S V]=svd(K_train);   LD_train=U*S^0.5;
% temp=cumsum(svd(LD_train));
% temp= temp./temp(end);
% dim=min(sum( (temp(2:end)- temp(1:end-1)) >0.001)+1, ...
%     max([3,sum(temp<dim_thresh)+1, sum(svd(LD_train)>0.1)]));
% LD_train=LD_train(:,1:dim);
% return;

%% calculate the dimensionality of the kernel matrix
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT:
% K_train       the kernel/grammian matrix of the training samples.  
% dim_thresh    the threshold on sigularvalue spectrum [0 1].
% OUTPUT:
% LD_train      the data represetatives of the training samples of the
%               learned dimensionality.
% dim           the learned dimensionality of the output manifold space.
function [LD_train, dim] = CalDim(K_train, dim_thresh)
[U S V]=svd(K_train);   LD_train=U*(S.^0.5);
sv = diag(S).^0.5;
temp=cumsum(sv);
sv = sv./temp(end);
temp= temp./temp(end);
dim=sum((temp<dim_thresh) | (sv >(1- dim_thresh)) );
LD_train=LD_train(:,1:dim);
return;