% Function that train and test the SDM manifold. The dimensionality of the
% output space is validated.
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
%       SDMStruct: the structure where all the training parameter and
%       trained data representatives are stored. The variables and
%       parameters are in the same name as in paper "The Support Vector Decomposition Machine"
%               + theta: sign(Z*theta) is the predicted label 
%               + Z, W: W*Z is used to approximated the data matrix 
%               + dim: the dimensionality of the projected space
%               + D: the weight of the hinge loss term 
%               + mu: the target classification margin.
%               + MisClassifiedIdx: the index of the validation set samples
%                       who are missed classified.
%               + LD_test: the low dimensional data representatives of the
%                       test set.
%               + Label: the predicted label annotation of the test set.
% TODO: the other output parameters are not really used. Need clean up 
function [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data,varargin)
num_train_fold=[];
if size(varargin)>0
    num_train_fold=varargin{1};
end
validation_fold=[];test_fold=[];SDMStruct=[];error_rate_SDM=[];
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
stc=struct('Z',[],'W',[],'theta',[],...
    'dim',10,'MisClassifiedIdx',[],...
    'LD_test',[],'Label',[],'D',[], 'mu', []);
SDMStruct=repmat(stc,size(exp_pair,1),length(num_fold));
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])
% load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_SDM_experiment.mat'])
if strcmp(name_data, 'face_grid')
    temp_data =norm_data;
    norm_data = zeros(size(norm_data,1)/4, size(norm_data,2));
    for i =1:size(norm_data,2)
        temp = temp_data(:,i);
        temp = imresize(reshape(temp,200,200), 0.5);
        norm_data(:,i)= temp(:);
    end
end
for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
         display(['================ #' num2str(i) ' fold ================'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_train_fold);
        %% SDM %%
        mu =1; % set as in the paper
        d_max =min( rank(train), 40);
        d_step = (d_max-1)/8;
        d_step = max(1, d_step);
        for d=1:8
            dim=round(d_step*d+1);
            if dim > size(train,2)
                continue;
            end
            parfor jj = 1:8
                D = 2^(2*(jj-1));
                [W, Z, theta, label_train,label_valid]=SDM(train, valid,group_train, dim, D, mu);
                W_valid{d,jj}=W;
                Z_valid{d,jj}=Z;
                theta_valid{d,jj}=theta;
                dim_valid{d,jj}=dim;
                D_valid{d,jj} = D;
                mu_valid{d,jj} = mu;
                error(d,jj)= sum(label_valid~=group_valid);
            end
            temp = error(d,:);
            temp(temp == 0) =max(temp); % do not choose the ones with 0 mis-detection, refer paper please.
            [error_temp(d) ix_temp(d)] =min(temp);
        end
        [~, d]=min(error_temp);
        SDMStruct(k,i).dim=dim_valid{d, ix_temp(d)}; SDMStruct(k,i).Z=Z_valid{d, ix_temp(d)};
        SDMStruct(k,i).W=W_valid{d, ix_temp(d)}; SDMStruct(k,i).theta=theta_valid{d, ix_temp(d)};
        SDMStruct(k,i).D=D_valid{d, ix_temp(d)}; SDMStruct(k,i).mu=mu_valid{d, ix_temp(d)};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   test               %%%%%%%%%%%%%%%%%%%
        W= SDMStruct(k,i).W;
        Z= SDMStruct(k,i).Z;
        theta= SDMStruct(k,i).theta;
        Z_test=(W'\test')';
        temp=Z_test*theta;
        temp(temp<0)=0;
        [temp ix]=sort(temp,2,'descend');
        label_test=zeros(size(Z_test,1),1);
        label_test(ix(:,1)==1)=1;
        label_test(ix(:,1)==2)=-1;
        %         [W, Z, theta, label_train,label_test]=SDM(train, test,group_train, dim);
        %         Z_res{k,i}=Z;
        %         W_res{k,i}=W;
        %         theta_res{k,i}=theta;
        error_rate_SDM(k,i)= sum(label_test~=group_test);
        SDMStruct(k,i).MisClassifiedIdx=test_idx(class_test~=group_test,:);
        save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) ...
            'fold_SDM_' num2str(num_train_fold) 'testfold_experiment.mat'],...
            'SDMStruct','error_rate_SDM', '-v7.3')
    end
end
return;

