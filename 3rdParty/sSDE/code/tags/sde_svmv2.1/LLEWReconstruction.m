% Compute and compare the Reconstruction error with local linear out of 
% sample extension method.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       name_data: data set name string
%       method: the method name of out of sample extension.
%       mode: out of sample extension mode name.
% OUTPUT
%       ReconError: the reconstruction error struct that contains 
%             the reconstruction error of the different dalgorithms
function [ReconError] = LLEWReconstruction(name_data, method, mode)
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
if strcmp(name_data, 'face_grid')
    temp_data =norm_data;
    norm_data = zeros(size(norm_data,1)/4, size(norm_data,2));
    for i =1:size(norm_data,2)
        temp = temp_data(:,i);
        temp = imresize(reshape(temp,200,200), 0.5);
        norm_data(:,i)= temp(:);
    end
end

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
            ReconError.SDESVM(i,:) = lleReconError_K(SDESVMStruct(i).K_train, ...
                train, SDESVMStruct(i).NN, SDESVMStruct(i).dim);
            ReconError.SDESVMKMF(i,:) = lleReconError_K(SDESVMKMFStruct(i).K_train, ...
                train, SDESVMKMFStruct(i).NN, SDESVMKMFStruct(i).dim);
            % SDE
            ReconError.SDE(i,:) = lleReconError_K(SDEStruct1(i).K_train, train, SDEStruct1(i).NN, SDEStruct1(i).dim);
            ReconError.SDEd(i,:) = lleReconError_K(SDEStruct1(i).K_train, train, SDEStruct2(i).NN, SDEStruct2(i).dim);
            
            % sLLE
            LD = slle(train', class_train, SLLEStruct1(i).NN, ...
                SLLEStruct1(i).dim, SLLEStruct1(i).alpha);
            ReconError.SLLE(i,:) = lleReconError(LD, train, SLLEStruct1(i).NN);
            LD = slle(train', class_train, SLLEStruct2(i).NN, ...
                SLLEStruct2(i).dim, SLLEStruct2(i).alpha);
            ReconError.SLLEd(i,:) = lleReconError(LD, train, SLLEStruct2(i).NN);
            
            % SDM
%             err = SDMStruct1(i).Z * SDMStruct1(i).W - train; 
%             err = sum(err.^2, 2)./size(train,1);%./ sum(train.^2, 2);
%             ReconError.SDM(i,:) = err;
%             err = SDMStruct2(i).Z * SDMStruct2(i).W - train; 
%             err = sum(err.^2, 2)./size(train,1);%./ sum(train.^2, 2); 
%             ReconError.SDMd(i,:) = err;
            ReconError.SDM(i,:) =lleReconError(SDMStruct1(i).Z, train, SDESVMKMFStruct(i).NN);
            ReconError.SDMd(i,:) =lleReconError(SDMStruct2(i).Z, train, SDESVMKMFStruct(i).NN);
        end
        sss =1;
        str =[name_data, ' & {$ ' num2str((mean(ReconError.SDE(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SDEd(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SDM(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SDMd(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SLLE(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SLLEd(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SDESVM(:))*sss),'%.4f') '$} '...
            ' & {$ ' num2str((mean(ReconError.SDESVMKMF(:))*sss),'%.4f') '$} \\'];
end  
display(str);
return;

function [err] = lleReconError_K(K, HD, NN, dim)
[U, D, V]=svd(K);
LD =U *(D.^0.5);
LD = LD(:,1:dim);
[err] = lleReconError(LD, HD, NN);
return;

function [err] = lleReconError(LD, HD, NN)
W=llew(LD', NN, inf);
W =W';
% err = norm(HD - W*HD, 'fro')^2;
err = HD - W*HD; 
err = sum(err.^2, 2);%./ sum(HD.^2, 2); ./size(HD,1)
return;
