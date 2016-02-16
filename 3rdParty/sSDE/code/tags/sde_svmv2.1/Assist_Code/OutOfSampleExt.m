%% out-of-sample extension
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUTs:
% K_train   the kernel/grammian matrix of the training samples. 
% St        the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% Sg        the high dimension test data, which is Ng-by-D, where Ng is the number
%           of samples
% y         the label of each sample {-1 +1}.
% NN        the number of local neightbors.
% dim_thresh    the threshold on sigularvalue spectrum [0 1].
% nCenters  the number of clusters,parameter for learning RBF mapping.
% LLM       the local linear projection matrices, for the RBF mapping.
% kmeans_space  the centers of each cluster, for RBF mapping.
% method    'RBF', NOT USED.
%           'MF'- learning the projection matrix for the testing samples as 
%                 in sSDE-l in the input space. Then compute the testing
%                 sample representatives in the output space with the
%                 learned projection matrix.
%           'LLW'-simply learn the local linear weight of each testing
%                 sample and apply these weights to compute the representatives 
%                 in the output low dimensional space. This is the method
%                 used in the ICML submission.
% mode      the embedding mode. 'SS' indicates the semisupervised mode,
%           which is not fully tested.
% OUTPUT:
% Training  the training data representatives in the output space, each
%           column is a embedded data sample.
% LD_test   the testing data representatives in the output space, each
%           column is a embedded data sample.
% cluster_label     the label of cluster for each sample
% LLM       the structure that store the local cluster information.
function [Training, LD_test, cluster_label,LLM] = ...
    OutOfSampleExt(K_train, St, Sg, y, NN,...
    dim_thresh, dim, nCenters,LLM, ...
    kmeans_space, method, mode)
LD_test=[]; cluster_label=[];
if sum(isnan(K_train(:)))>0
    Training =[];
    LD_test =[];
    cluster_label =[];
    LLM =[];
    return;
end
% calculate dimensionality
if isempty(dim)
    [LD, dim] = CalDim(K_train, dim_thresh);
    dim
else
    [U S V]=svd(K_train);   LD_train=U*S^0.5;
    LD=LD_train(:,1:dim);
end
Nt=size(St,1);
Training=LD(1:Nt,1:dim);

% embedding the out of sample test data and classify
if strcmp(mode, 'SS')
    LD_test = LD(Nt+1:end,1:dim);
else    
    if strcmp(method, 'RBF')
        [LD_test,cluster_label,LLM] = ...
            NLMap(St,Training,Sg, y, nCenters, LLM, kmeans_space);
    elseif strcmp(method, 'MF')
        [LD_test] = ISDEMap(St,Training,Sg, NN);
    elseif strcmp(method, 'LLW')
        [LD_test] = ISDEMap(St,Training,Sg, NN, method);
    end
end
return;