% Parameter setting for different dataset. Since the incomplete
% parameter cross-validation is applied for sSDE algorithm, the initial
% parameters is set up here.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       name_data: the string for data set name
% OUTPUT
%       C: the weight on error term.
%       lam: the weight between trace(K) and the sumation of margin and
%           errors
%       NN: the size of local neighborhood.
%       eps: the relaxation ratio of the local isometric constraints
%       nCenters(NOT USED): number of clusters for RBF extension
function [C, lam, NN, eps, nCenters]=InitPar(name_data)
C=16; NN=4;eps=0.1; lam=1e-3; nCenters=20;
if strcmp(name_data, 'mug_silhouette')
    C=16; NN=6;eps=0.1; lam=1e-3; nCenters=20;
elseif strcmp(name_data, 'mouse_silhouette')
    C=16; NN=6;eps=0.1; lam=1e-3; nCenters=20;
elseif strcmp(name_data, 'stapler_silhouette')
    C=16; NN=6;eps=0.1; lam=1e-3; nCenters=20;
% elseif strcmp(name_data, 'PAL')
%     C=4; NN=4;eps=0.1; lam=1e-3; nCenters=20;
elseif strcmp(name_data, 'Ionosphere')
    C=16; NN=6;eps=0.1; lam=1e-3; nCenters=20;
% elseif strcmp(name_data, 'sonar')
%     C=16; NN=5;eps=0.1; lam=1e-3; nCenters=20;
end