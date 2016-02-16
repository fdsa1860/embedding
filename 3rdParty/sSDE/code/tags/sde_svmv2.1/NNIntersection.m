% Compute and display the intersection of the K nearest neiborhood set between 
% input and output space.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       name_data: data set name string
%       method: the method name of out of sample extension.
%       mode: out of sample extension mode name.
% OUTPUT
%       ReconError: the ratio of the intersectionof the K nearest neiborhood 
%           set between input and output space.
function [ReconError] = NNIntersection(name_data, method, mode)
num_test_fold =[];
class =[];
load([name_data '_nomalized.mat'])
protocol;
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])
SDESVMStruct =[];
SDESVMKMFStruct =[];
SDEStruct1 = [];
SDEStruct2 = [];
SDMStruct1 = [];
SDMStruct2 =[];
SLLEStruct1 = [];
SLLEStruct2 =[];


Nl=round(0.3*num_samp*4);
if strcmp(name_data, 'PAL')
    Nl=round(0.3*num_samp*4*2);
end
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDE_'  num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
    'SDEStruct');
SDEStruct1 = SDEStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDEv2_' num2str(Nl) '_' method '_' mode  '_experiment.mat'],...
    'SDEStruct');
SDEStruct2 = SDEStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDM_'  num2str(num_test_fold) 'testfold_experiment.mat'],...
    'SDMStruct');
SDMStruct1 = SDMStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDMv2_' num2str(Nl) '_experiment.mat'],...
    'SDMStruct');
SDMStruct2 = SDMStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_sLLEv1_experiment.mat'],...
    'SLLEStruct');
SLLEStruct1 = SLLEStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_sLLEv2_' num2str(Nl) ' _experiment.mat'],...
    'SLLEStruct');
SLLEStruct2 = SLLEStruct;
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDESVM_'  num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
    'SDESVMStruct');
%----------------------------------------------------------
load([name_data '_' num2str(num_samp) 'sample' '_' ...
    num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
    method '_' mode '_experiment.mat'],...
    'SDESVMKMFStruct');
if strcmp(name_data, 'PAL')
    num_samp = num_samp * 2;
end
for k=1:size(exp_pair,1);
        pos_class=exp_pair(k,1);
        neg_class=exp_pair(k,2);
        for i=1:num_fold
            [train, test, valid, group_train, group_test,group_valid, ...
                class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
                FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
            % SDESVM
            ReconError.SDESVM(i,:) = NNIntersec_K(SDESVMStruct(i).K_train, ...
                train, SDESVMStruct(i).NN, SDESVMStruct(i).dim);
            ReconError.SDESVMKMF(i,:) = NNIntersec_K(SDESVMKMFStruct(i).K_train, ...
                train, SDESVMKMFStruct(i).NN, SDESVMKMFStruct(i).dim);
            % SDE
            ReconError.SDE(i,:) = NNIntersec_K(SDEStruct1(i).K_train, train, SDEStruct1(i).NN, SDEStruct1(i).dim);
            ReconError.SDEd(i,:) = NNIntersec_K(SDEStruct2(i).K_train, train, SDEStruct2(i).NN, SDEStruct2(i).dim);
            
            % sLLE
            LD = slle(train', class_train, SLLEStruct1(i).NN, ...
                SLLEStruct1(i).dim, SLLEStruct1(i).alpha);
            ReconError.SLLE(i,:) = NNIntersec(LD, train, SLLEStruct1(i).NN);
            LD = slle(train', class_train, SLLEStruct2(i).NN, ...
                SLLEStruct2(i).dim, SLLEStruct2(i).alpha);
            ReconError.SLLEd(i,:) = NNIntersec(LD, train, SLLEStruct2(i).NN);
            
            % SDM
            ReconError.SDM(i,:) = NNIntersec(SDMStruct1(i).Z, train, SDESVMKMFStruct(i).NN);
            ReconError.SDMd(i,:) = NNIntersec(SDMStruct2(i).Z, train, SDESVMKMFStruct(i).NN);
        end
        str =[name_data, ' & {$ ' num2str(round(mean(ReconError.SDE(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SDEd(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SDM(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SDMd(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SLLE(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SLLEd(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SDESVM(:))*100)/100) '$} '...
            ' & {$ ' num2str(round(mean(ReconError.SDESVMKMF(:))*100)/100) '$} \\'];
end
display(str);
return;

function [err] = NNIntersec_K(K, HD, NN, dim)
[U, D, V]=svd(K);
LD =U *(D.^0.5);
LD = LD(:,1:dim);
[err] = NNIntersec(LD, HD, NN);
return;

function [err] = NNIntersec(LD, HD, NN)
[neighbors_LD]= NNeighbour(LD', NN);
[neighbors_HD]= NNeighbour(HD', NN);
parfor i =1:size(LD,1)
    c = intersect(neighbors_LD(:,i)', neighbors_HD(:,i)');
    err(i,1) = 1-length(c)/NN;
end
return;

function [neighbors]= NNeighbour(X, K)
[D,N] = size(X);
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
% NEIGHBORS
for i=1:N
    distance(i,i)=max(distance(i,:));
end
[sorted,index] = sort(distance);
neighbors = index(1:K,:);
return;