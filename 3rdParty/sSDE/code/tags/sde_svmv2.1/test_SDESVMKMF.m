% Function to cross-validate and test the sSDE_l.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
%       name_data: data set name string
%       dim_thresh: the threshold to determine the dimensionality of the
%                   embedded manifold.
%       method: the method name of out of sample extension.
%       mode: out of sample extension mode name.
%   
% OUTPUT
%       SDESVMKMFStruct: the structure where all the training parameter and
%       trained data representatives are stored.
%               + K_train: the Grammian matrix of the training data
%                       samples.
%               + NN: the cross-valided number of neighborhood
%               + Nl: the number of landmarks
%               + lambda: the weight ratio parameter [0, 1] for the trace(K).
%               + C: the cross-valided weight of the SVM error term
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
%               + t: the supreme limit of the SVM objective 
%               + cvx_result: the optimum value of the cvx.
%               + run_time: the time cvx solve the problem.
% TODO: the other output parameters are not really used. Need clean up      
function [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVMKMF] ...
    = test_SDESVMKMF(name_data,dim_thresh, method, mode , varargin)
display(['============'  name_data '=============='])
display('============SDESVMKMFM==============')
num_test_fold=[];
if size(varargin)>0
    num_test_fold=varargin{1};
end
validation_fold=[];test_fold=[];SDESVMKMFStruct=[];error_rate_SDESVMKMF=[];
norm_data=[];
datset=[];
test_fold=[];
validation_fold=[];
class=[];
num_samp=20;
load([name_data '_nomalized.mat']) % read the data
protocol;
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])
stc=struct('u',[],'gamma',[],'lambda',[],'C',[],'NN',[],'eps',[],...
    'dim',10,'K_train',[],'nCenters',[],'MisClassifiedIdx',[]...
    ,'LLM',[],'LD_test',[],'Label',[],'indl',[] ,'t',[]...
    ,'v1',[] ,'v2',[] ,'v3',[] ,'cvx_result',[]);
SDESVMKMFStruct=repmat(stc,size(exp_pair,1),num_fold);
fname=[];
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDESVM_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end
% fname=[name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDE_experiment.mat'];
% if exist(fname)==2
%     load(fname);
% end

% fname_temp=[name_data '_' num2str(num_samp) 'sample' '_' ...
%     num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
%     method '_' mode '_experiment.mat'];
% if exist(fname_temp)==2
%     load(fname_temp);
% end
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        if i<=size(SDESVMKMFStruct,2) && ~isempty(SDESVMKMFStruct(k,i).K_train)
            continue;
        end
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        error_rate_temp=[];
        %% SDESVM %%
        Nl=round(0.3*size(train,1));
        [C, lam, NN, eps, nCenters]=InitPar(name_data);
        SDESVMKMFStruct(k,i).eps=eps;
        SDESVMKMFStruct(k,i).Nl = Nl;
        %% valid lam
        display('----------------- valid lam ---------------------')
        error_rate_temp=[];
        parfor j=1:8
            lam=1/(1+10^(j-4));
            [K_train eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdual(train, valid,NN,Nl, eps,group_train,lam, C, mode);
            K_valid1{j,1}=K_train; t_valid{j,1}=t;
            lam_valid{j,1} =lam;
        end
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
            'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl','K_valid1', '-v7.3')
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
        SDESVMKMFStruct(k,i).lambda=lam;
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
            'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl','K_valid1', '-v7.3')
                %% valid C
        display('----------------- valid C ---------------------')
        error_rate_temp=[];
        lam = SDESVMKMFStruct(k,i).lambda;
        parfor c=1:8
            C=2^(2*(c-1));
            [K_train eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdual(train, valid,NN,Nl, eps,group_train,lam, C, mode);
            K_valid2{c,1}=K_train;  
            C_valid{c, 1} =C;
        end
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
            'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl','K_valid2', '-v7.3')
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
        SDESVMKMFStruct(k,i).C=C;
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
            'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl','K_valid2', '-v7.3')
    end
end
%% valid NN
display('----------------- valid NN Embedding Stage---------------------')
time_valid3 =zeros(4, num_fold);
NN_valid =zeros(4, num_fold);
K_valid3 =cell(4, num_fold);
indl_valid3=cell(4, num_fold);
t_valid3 =cell(4, num_fold);
v1_valid3 =cell(4, num_fold);
v2_valid3 =cell(4, num_fold);
v3_valid3 =cell(4, num_fold);
cvx_result_valid3 = cell(4, num_fold);
%%
for j=1:4
    NN=3+j;
    for k=1:size(exp_pair,1);
        pos_class=exp_pair(k,1);
        neg_class=exp_pair(k,2);
        display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
        parfor i=1:num_fold
            display(['================  #' num2str(i) 'fold ================'])
            if i<=size(SDESVMKMFStruct,2) && ~isempty(SDESVMKMFStruct(k,i).K_train)
                continue;
            end
            [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
                FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
            
            time1 = cputime;
            [K_train eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] ...
                = SDESVMKMFdual(train, valid,NN,Nl, SDESVMKMFStruct(k,i).eps,group_train,...
                SDESVMKMFStruct(k,i).lambda, SDESVMKMFStruct(k,i).C, mode);
            time2 = cputime;
            time_valid3(j,i) = time2 -time1;
            
            K_valid3{j,i}=K_train;
            NN_valid(j,i) = NN;
            indl_valid3{j,i}=indl;
            t_valid3{j,i} =t;
            v1_valid3{j,i} =v1;
            v2_valid3{j,i} =v2;
            v3_valid3{j,i} =v3;
            cvx_result_valid3{j,i} = cvx_result;            
        end
    end
    
    save([name_data '_' num2str(num_samp) 'sample' '_' ...
        num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
        method '_' mode '_experiment.mat'],...
        'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl',...
        'K_valid3','NN_valid','t_valid3','indl_valid3',...
        'time_valid3','v1_valid3','v2_valid3','v3_valid3','cvx_result_valid3', '-v7.3')

end
%%
u_valid3 =cell(size(K_valid3, 1), 1);
gamma_valid3 =cell(size(K_valid3, 1), 1);
dim_valid3 =cell(size(K_valid3, 1), 1);
display('----------------- valid NN Validation Stage---------------------')
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        if i<=size(SDESVMKMFStruct,2) && ~isempty(SDESVMKMFStruct(k,i).K_train)
            continue;
        end
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        eps = SDESVMKMFStruct(k,i).eps;
        C = SDESVMKMFStruct(k,i).C;
        lam = SDESVMKMFStruct(k,i).lambda;
        NN  = NN_valid(j,i);
        nCenters=min(floor(size(train,1)/2),20);

        error_rate_temp= zeros(size(K_valid3, 1),1);
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
        SDESVMKMFStruct(k,i).NN=NN;
        SDESVMKMFStruct(k,i).K_train=K_valid3{j,i};    
        SDESVMKMFStruct(k,i).u=u_valid3{j,1};  
        SDESVMKMFStruct(k,i).gamma=gamma_valid3{j,1};
        SDESVMKMFStruct(k,i).dim=dim_valid3{j,1};
        SDESVMKMFStruct(k,i).indl = indl_valid3{j,i};
        SDESVMKMFStruct(k,i).v1 = v1_valid3{j,i};
        SDESVMKMFStruct(k,i).v2 = v2_valid3{j,i};
        SDESVMKMFStruct(k,i).v3 = v3_valid3{j,i};
        SDESVMKMFStruct(k,i).t = t_valid3{j,i};
        SDESVMKMFStruct(k,i).cvx_result = cvx_result_valid3{j,i};
        SDESVMKMFStruct(k,i).run_time=time_valid3(j,i);
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
        %             u=SDESVMKMFStruct(k,i).u;  gamma=SDESVMKMFStruct(k,i).gamma;
        %             label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
        %             LLM_valid{n,1}=LLM_temp;
        %             error_rate_temp(n,1)= sum(label_valid~=group_valid);
        %             nCenters_valid{n, 1} =nCenters;
        %         end
        %         display([ 'error rate: ' num2str(error_rate_temp)])
        %         [~,n]=min(error_rate_temp); nCenters=nCenters_valid{n, 1};
        %         SDESVMKMFStruct(k,i).nCenters=nCenters;
        %         SDESVMKMFStruct(k,i).LLM=LLM_valid{n,1};
        
        save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
            'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl',...
            'K_valid3','NN_valid','t_valid3','indl_valid3',...
            'time_valid3','v1_valid3','v2_valid3','v3_valid3','cvx_result_valid3', '-v7.3')
    end
end


for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%         test         %%%%%%%%%%%%%%%%%%%
    parfor i=1:num_fold
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        K_train=SDESVMKMFStruct(k,i).K_train;
        NN=SDESVMKMFStruct(k,i).NN;
        dim=SDESVMKMFStruct(k,i).dim;
        LLM=SDESVMKMFStruct(k,i).LLM;
        C=SDESVMKMFStruct(k,i).C;
        if strcmp(mode, 'SS')
            lam=SDESVMKMFStruct(k,i).lambda;
            [K_train eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdual(train, test,NN,Nl, eps,group_train,lam, C, mode);
        end
        [LD_train, LD_test] = OutOfSampleExt(K_train, train, test, group_train, NN,...
            [], dim,[], LLM, 'LD', method, mode);
        [u, gamma, e]=train_svm(LD_train, group_train, C);
        label_test =test_svm(LD_train, group_train, LD_test, u, gamma);
        error_rate_SDESVMKMF(k,i)= sum(label_test~=group_test);
        SDESVMKMFStruct(k,i).MisClassifiedIdx=test_idx(label_test~=group_test,:);
        SDESVMKMFStruct(k,i).LD_test=LD_test;
        SDESVMKMFStruct(k,i).Label=label_test;
    end
end
[sum(error_rate_SDESVMKMF,2)]
save([name_data '_' num2str(num_samp) 'sample' '_' ...
            num2str(num_fold) 'fold_SDESVMKMF_' num2str(Nl) '_' num2str(num_test_fold) 'testfold_'...
            method '_' mode '_experiment.mat'],...
    'SDESVMKMFStruct','error_rate_SDESVMKMF','dim_thresh','Nl', '-v7.3')
%%

return;
function [K eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdual(St,Sg, NN,Nl,eps,y,lam1, lam2, mode,varargin)
% [K eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdualv0(St,Sg, NN,Nl,eps,y,lam1, lam2, mode,varargin);
[K eta indl cvx_optval cvx_status t v1 v2 v3 cvx_result] = SDESVMKMFdualv1(St,Sg, NN,Nl,eps,y,lam1, lam2, mode,varargin);
return; 

