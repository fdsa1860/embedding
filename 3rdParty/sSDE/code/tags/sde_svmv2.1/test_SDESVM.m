% Function to cross-validate and test the sSDE.
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
    = test_SDESVM(name_data,dim_thresh, method, mode , varargin)
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
stc=struct('u',[],'gamma',[],'lambda',[],'C',[],'NN',[],'eps',[],...
    'dim',[],'K_train',[],'nCenters',[],'MisClassifiedIdx',[]...
    ,'LLM',[],'LD_test',[],'Label',[],'run_time',[]);
SDESVMStruct=repmat(stc,size(exp_pair,1),num_fold);
fname=[];
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDESVM_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDE_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        if i<=size(SDESVMStruct,2) && ~isempty(SDESVMStruct(k,i).K_train)
            continue;
        end
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        %% SDESVM %%
        C=16; NN=4;eps=0.1; lam=1e-3; nCenters=min(floor(size(train,1)/2),20);
        SDESVMStruct(k,i).eps=eps;
       
        %% valid lam
        if i<=size(SDESVMStruct,2) && isempty(SDESVMStruct(k,i).lambda)
            display('----------------- valid lam ---------------------')
            error_rate_temp=[];
            parfor j=1:8
                pause(10*(8-j+1));
                lam=1/(1+10^(j-4));
                [K_train eta cvx_optval cvx_status t] = SDESVMdualv0(train, valid,NN,eps,group_train,lam, C, mode);
                K_valid1{j,1}=K_train; t_valid{j,1}=t;
                lam_valid{j,1} =lam;
            end
            save([name_data '_' num2str(num_samp) 'sample' '_' ...
                num2str(num_fold) 'fold_SDESVM_' num2str(num_test_fold) 'testfold_'...
                method '_' mode '_experiment.mat'], 'SDESVMStruct','K_valid1','-v7.3')
            parfor j=1:8
                K_train=K_valid1{j,1};
                [LD_train, LD_valid] = OutOfSampleExt(K_train, train, valid, group_train,NN,...
                    dim_thresh,[], nCenters, [], 'LD', method, mode);
                if isempty(LD_train)
                    error_rate_temp(j,1) = size(valid,1);
                    continue;
                end
                [u, gamma, e]=train_svm(LD_train, group_train, C);
                label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
                error_rate_temp(j,1)= sum(label_valid~=group_valid);
            end
            error_rate_temp(error_rate_temp~=min(error_rate_temp))=inf;
            [~, j]=min(error_rate_temp);    lam=lam_valid{j,1};
            SDESVMStruct(k,i).lambda=lam;
        end
        %% valid C
        if i<=size(SDESVMStruct,2) && isempty(SDESVMStruct(k,i).C)
            display('----------------- valid C ---------------------')
            error_rate_temp=[];
            lam = SDESVMStruct(k,i).lambda;
            parfor c=1:8
                pause(10*(8-c+1));
                C=2^(2*(c-1));
                [K_train eta cvx_optval cvx_status] = SDESVMdualv0(train, valid,NN,eps,group_train,lam, C, mode);
                K_valid2{c,1}=K_train;
                C_valid{c, 1} =C;
            end
            save([name_data '_' num2str(num_samp) 'sample' '_' ...
                num2str(num_fold) 'fold_SDESVM_' num2str(num_test_fold) 'testfold_'...
                method '_' mode '_experiment.mat'], 'SDESVMStruct','K_valid2','-v7.3')
            parfor c=1:8
                K_train=K_valid2{c,1};
                C=C_valid{c, 1};
                [LD_train, LD_valid] = OutOfSampleExt(K_train, train, valid, group_train,NN,...
                    dim_thresh,[], nCenters, [], 'LD', method, mode);
                if isempty(LD_train)
                    error_rate_temp(c,1) = size(valid,1);
                    continue;
                end
                [u, gamma, e]=train_svm(LD_train, group_train, C);
                label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
                error_rate_temp(c,1)= sum(label_valid~=group_valid);
            end
            error_rate_temp(error_rate_temp~=min(error_rate_temp))=inf;
            [~, c]=min(error_rate_temp);    C=C_valid{c, 1};
            SDESVMStruct(k,i).C=C;
        end
    end
end
%% valid NN
display('----------------- valid NN Embedding Stage---------------------')
t_valid3 =zeros(4, num_fold);
NN_valid =zeros(4, num_fold);
K_valid3 =cell(4, num_fold);
%%
for j=1:4
    NN=3+j;
    for k=1:size(exp_pair,1);
        pos_class=exp_pair(k,1);
        neg_class=exp_pair(k,2);
        display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
        parfor i=1:num_fold
            display(['================  #' num2str(i) 'fold ================'])
            if i<=size(SDESVMStruct,2) && ~isempty(SDESVMStruct(k,i).K_train)
                continue;
            end
            [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
                FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
            
            time1 = cputime;
            
            [K_train eta cvx_optval cvx_status] = SDESVMdualv0(train, valid,NN,...
                eps,group_train,SDESVMStruct(k,i).lambda, SDESVMStruct(k,i).C, mode);
            
            time2 = cputime;
            t_valid3(j,i) = time2 -time1;
            
            K_valid3{j,i}=K_train;
            NN_valid(j,i) = NN;
        end
    end
    
    save([name_data '_' num2str(num_samp) 'sample' '_' ...
        num2str(num_fold) 'fold_SDESVM_' num2str(num_test_fold) 'testfold_'...
        method '_' mode '_experiment.mat'], 'SDESVMStruct','K_valid3','t_valid3','NN_valid','-v7.3')
end
%%
display('----------------- valid NN Validation Stage---------------------')
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        if i<=size(SDESVMStruct,2) && ~isempty(SDESVMStruct(k,i).K_train)
            continue;
        end
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        eps = SDESVMStruct(k,i).eps;
        C = SDESVMStruct(k,i).C;
        lam = SDESVMStruct(k,i).lambda;
        NN  = NN_valid(j,i);
        nCenters=min(floor(size(train,1)/2),20);

        error_rate_temp=[];
        parfor j=1:size(K_valid3, 1)
            K_train=K_valid3{j,i};
            [LD_train, LD_valid] = OutOfSampleExt(K_train, train, valid, group_train,NN,...
                dim_thresh,[], nCenters, [], 'LD', method, mode);
            if isempty(LD_train)
                error_rate_temp(j,1) = size(valid,1);
                continue;
            end
            [u, gamma, e]=train_svm(LD_train, group_train, C);
            u_valid3{j,1}=u; gamma_valid3{j,1}=gamma; dim_valid3{j,1}=size(LD_train,2);
            label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
            error_rate_temp(j,1)= sum(label_valid~=group_valid);
        end
        display(num2str(error_rate_temp));
        error_rate_temp(error_rate_temp~=min(error_rate_temp))=inf;
        [~, j]=min(error_rate_temp);    NN=NN_valid(j,i);
        SDESVMStruct(k,i).NN=NN;
        SDESVMStruct(k,i).K_train=K_valid3{j,i};    
        SDESVMStruct(k,i).u=u_valid3{j,1};  
        SDESVMStruct(k,i).gamma=gamma_valid3{j,1};
        SDESVMStruct(k,i).dim=dim_valid3{j,1};
        SDESVMStruct(k,i).run_time=t_valid3(j,i);
%         %% valid nCenters
%         error_rate_temp=[];
%         K_train=K_valid3{j,1};
%         parfor n=1:10
%             nCenters=3*n;
%             if nCenters>=size(train,1)/2
%                 continue;
%             end
%             [LD_train, LD_valid, ~, LLM_temp] = OutOfSampleExt(K_train, train, valid, group_train,NN,...
%                 dim_thresh,[], nCenters, [], 'LD', method, mode);
%             if isempty(LD_train)
%                 error_rate_temp(n,1) = size(valid,1);
%                 continue;
%             end
%             u=SDESVMStruct(k,i).u;  gamma=SDESVMStruct(k,i).gamma;
%             label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
%             LLM_valid{n,1}=LLM_temp;
%             error_rate_temp(n,1)= sum(label_valid~=group_valid);
%             nCenters_valid{n, 1} =nCenters;
%         end
%         display([ 'error rate: ' num2str(error_rate_temp)])
%         [~,n]=min(error_rate_temp); nCenters=nCenters_valid{n, 1};
%         SDESVMStruct(k,i).nCenters=nCenters;
%         SDESVMStruct(k,i).LLM=LLM_valid{n,1};
       
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVM_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'], 'SDESVMStruct','t_valid3','-v7.3')
    end
end

%%
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%         test         %%%%%%%%%%%%%%%%%%%
    parfor i=1:num_fold
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        K_train=SDESVMStruct(k,i).K_train;
        NN=SDESVMStruct(k,i).NN;
        dim=SDESVMStruct(k,i).dim;
        LLM=SDESVMStruct(k,i).LLM;
        C=SDESVMStruct(k,i).C;
        if strcmp(mode, 'SS')
            lam=SDESVMStruct(k,i).lambda;
            [K_train eta cvx_optval cvx_status] = SDESVMdualv0(train, test,NN,eps,group_train,lam, C, mode);
        end
        [LD_train, LD_test] = OutOfSampleExt(K_train, train, test, group_train, NN,...
            [], dim,[], LLM, 'LD', method, mode);
        [u, gamma, e]=train_svm(LD_train, group_train, C);
        label_test =test_svm(LD_train, group_train, LD_test, u, gamma);
        error_rate_SDESVM(k,i)= sum(label_test~=group_test);
        SDESVMStruct(k,i).MisClassifiedIdx=test_idx(label_test~=group_test,:);
        SDESVMStruct(k,i).LD_test=LD_test;
        SDESVMStruct(k,i).Label=label_test;
    end
end
[sum(error_rate_SDESVM,2)]

save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDESVM_' num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
    'SDESVMStruct','error_rate_SDESVM','t_valid3', '-v7.3')
%%
return;


