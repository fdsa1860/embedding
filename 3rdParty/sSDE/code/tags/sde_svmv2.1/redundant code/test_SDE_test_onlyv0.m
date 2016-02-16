function [validation_fold,test_fold,SDEStruct,error_rate_SDE] = test_SDE_test_onlyv0(name_data,dim_thresh,method, mode , varargin)
num_test_fold=[];
if size(varargin)>0
    num_test_fold=varargin{1};
end
validation_fold=[];test_fold=[];SDEStruct=[];error_rate_SDE=[];
norm_data=[];
datset=[];
test_fold=[];
validation_fold=[];
class=[];
load([name_data '_nomalized.mat'])
protocol;
stc=struct('u',[],'gamma',[],'lambda',[],'C',[],'NN',[],'eps',[],...
    'dim',10,'K_train',[],'nCenters',[],'MisClassifiedIdx',[]...
    ,'LLM',[],'LD_test',[],'Label',[]);
SDEStruct=repmat(stc,size(exp_pair,1),length(num_fold));
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])
% load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDESVM_experiment.mat'])
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDE_'  num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat']);
%% SVM TRAINING
valid_sz =[4 7];
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        eps =0.1;
        nCenters =10;  % Not used here
        for c=1:valid_sz(2)
            C=2^(2*(c-1));
            parfor j=1:valid_sz(1)
                NN=4+j;
                [LD_train, LD_valid, ~, LLM_temp] = OutOfSampleExt(K_valid3{j,i}, train, valid, group_train, NN,...
                    dim_thresh, [],  nCenters,[], 'LD', method, mode);
                [u, gamma, e, score(j,c), margin(j,c)]=train_svm(LD_train, group_train, C);
                u_valid3{j,c}=u; gamma_valid3{j,c}=gamma;
                dim_valid3{j,c}=size(LD_train,2);
                LLM_valid3{j,c}=LLM_temp;
                label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
                error_rate_temp(j,c)= sum(label_valid~=group_valid);
            end
        end
        temp=zeros(size(error_rate_temp));   temp(error_rate_temp~=min(error_rate_temp(:)))=inf;
        ix=find(temp==min(temp(:)));
        [j c]=ind2sub(valid_sz,ix(1));
        NN=4+j;  C=2^(2*(c-1));
        SDEStruct(k,i).NN=NN; SDEStruct(k,i).C=C;
        SDEStruct(k,i).K_train=K_valid3{j,i}; 
        SDEStruct(k,i).u=u_valid3{j,c};  
        SDEStruct(k,i).gamma=gamma_valid3{j,c};
        SDEStruct(k,i).dim=dim_valid3{j,c};
        SDEStruct(k,i).LLM=LLM_valid3{j,c};
        SDEStruct(k,i).eps=eps;
        save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
            'fold_SDE_'  num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
            'SDEStruct','K_valid3','t_valid3', '-v7.3')
    end
    parfor i=1:num_fold
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_test_fold);
%         display(['================ test ================'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   test               %%%%%%%%%%%%%%%%%%%
        K_train=SDEStruct(k,i).K_train;
        NN=SDEStruct(k,i).NN;
        dim=SDEStruct(k,i).dim;
        C=SDEStruct(k,i).C;
        LLM=SDEStruct(k,i).LLM;
        
        [~, dim] = CalDim(K_train, dim_thresh); % correction for former bug 08-22-2013
        
        if strcmp(mode, 'SS')
            X=[train;test];
            [LD, ~] = SDEcvx(X, eps, NN);
            SDEStruct(k,i).K_test=LD*LD';
            K_train = SDEStruct(k,i).K_test;
        end
        
        [LD_train, LD_test] = OutOfSampleExt(K_train, train, test, group_train, NN,...
            dim_thresh, [],[], LLM, 'LD', method, mode);
        [u, gamma, e]=train_svm(LD_train, group_train, C);
        label_test =test_svm(LD_train, group_train, LD_test, u, gamma);
        error_rate_SDE(k,i)= sum(label_test~=group_test);
        SDEStruct(k,i).u = dim;
        SDEStruct(k,i).u = u;
        SDEStruct(k,i).gamma = gamma;
        SDEStruct(k,i).MisClassifiedIdx=test_idx(label_test~=group_test,:);
        SDEStruct(k,i).LD_test=LD_test;
        SDEStruct(k,i).Label=label_test;
    end
    Nl=round(0.3*size(SDEStruct(1,1).K_train,1));
    save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
    'fold_SDE_'  num2str(num_test_fold) 'testfold_' method '_' mode '_experiment.mat'],...
    'SDEStruct','error_rate_SDE','K_valid3','t_valid3', '-v7.3')
end
sum(error_rate_SDE,2)
return;
