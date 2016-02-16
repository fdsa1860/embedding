% Function to cross-validate and test the sSDE on only the landmark samples.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       name_data: data set name string
%       dim_thresh: 
%       method: the method name of out of sample extension.
%       mode: out of sample extension mode name.
%   
% OUTPUT
%       SDESVMKMFStruct: the structure where all the training parameter and
%       trained data representatives are stored.
%               + K_train: the Grammian matrix of the training data
%                       samples.
%               + NN: the cross-valided number of neighborhood
%               + lambda: the weight ratio parameter [0, 1] for the trace(K).
%               + C: the cross-valided the weight of the SVM error term
%               + u: the dual variable of the SVM
%               + gamma: the offset of the boundary
%               + dim: computed dimensionality of the manifold.
%               + LLM: local projection matrix
%               + eps: local geometric constraint relaxation ratio.
%               + MisClassifiedIdx: the index of the validation set samples
%                       who are missed classified.
%               + LD_test: the low dimensional data representatives of the
%                       test set.
%               + Label: the predicted label annotation of the test set.
%               + v1, v2, v3: the dual variable in paper.
%               + t: the supreme bond of the SVM objective 
%               + cvx_result: the optimum value of the cvx.
%               + run_time: the time cvx solve the problem.
% TODO: the other output parameters are not really used. Need clean up  
function [validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] ...
    = test_SDESVM_nl_only(name_data, fname, dim_thresh, method, mode , varargin)
num_test_fold=[];
if size(varargin)>0
    num_test_fold=varargin{1};
end
validation_fold=[];test_fold=[];SDESVMStruct=[];error_rate_SDESVM1=[];


norm_data=[];
datset=[];
test_fold=[];
validation_fold=[];
class=[];
num_samp=20;
% name_data='usps';
% name_data='PAL';
load([name_data '_nomalized.mat'])
protocol;
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])
load(fname);
% fname=[];
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDESVM_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDE_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
                'fold_SDESVMKMF_nl_' num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat']);
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        error_rate_temp=[];
        %% SDESVM %%
        nCenters =20; % not used
        eps = SDESVMStruct(k,i).eps;
        NN = SDESVMStruct(k,i).NN;
        C = SDESVMStruct(k,i).C;
        lam = SDESVMStruct(k,i).lambda;
        %% use the parameter from SDESVMStruct to test SDESVMKMF
        for m = 1:8
            Nl = floor(m*0.05*size(train,1));
            load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
                'fold_SDESVMKMF_nl_' num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat']);
            parfor j=1:8
                training = train(idx_landmarks{i,m,j},:);
                group_training = group_train(idx_landmarks{i,m,j},:);
                % C = C*size(train,1)/ length(idx_landmarks{i,m,j});
                % lam = lam*size(train,1)/ length(idx_landmarks{i,m,j});
                t1 = cputime;
                [K_train eta cvx_optval cvx_status t] = SDESVMdualv0(training, valid,NN,eps,group_training,lam, C, mode);
                t2 = cputime;
                run_time{i,m,j} = t2-t1;
                K{i,m,j}=K_train;
            end
            parfor j=1:8
                training = train(idx_landmarks{i,m,j},:);
                group_training = group_train(idx_landmarks{i,m,j},:);
                K_train=K{i,m,j};
                [LD_train, LD_test] = OutOfSampleExt(K_train, training, test, group_training, NN,...
                    dim_thresh,[], nCenters, [], 'LD', method, mode);
                [u, gamma, e]=train_svm(LD_train, group_training, C);
                label_test =test_svm(LD_train, group_training, LD_test, u, gamma);
                error_rate_SDESVM_nl(i,m,j)= sum(label_test~=group_test)/length(group_test);
            end
            save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
                'fold_SDESVM_nl_only_' num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
                'SDESVMStruct','error_rate_SDESVM_nl', 'K', 'run_time', '-v7.3')
            [mean(error_rate_SDESVM_nl,2)]
        end
    end
end

%%

return;


