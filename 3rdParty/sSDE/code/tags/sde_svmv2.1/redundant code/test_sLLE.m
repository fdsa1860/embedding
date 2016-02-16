function [validation_fold,test_fold,SLLEStruct,error_rate_sLLE] = test_sLLE(name_data, method, mode , varargin)
num_train_fold=[];
if size(varargin)>0
    num_train_fold=varargin{1};
end
validation_fold=[];test_fold=[];SLLEStruct=[];error_rate_sLLE=[];
SDESVMStruct=[];
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
SLLEStruct=repmat(struct('alpha', [],'NN',[],'dim',10,'C',[],'u',[],'gamma',[],'MisClassifiedIdx',[]),size(exp_pair,1),length(num_fold));
load([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold.mat'])


for k=1:size(exp_pair,1);
    pos_class=exp_pair(k,1);
    neg_class=exp_pair(k,2);
    display(['================ ' num2str(pos_class) ' vs ' num2str(neg_class) ' ================'])
    for i=1:num_fold
        display(['================  #' num2str(i) 'fold ================'])
        [train, test, valid, group_train, group_test,group_valid, class_train, class_test, class_valid, train_idx, test_idx, ~]= ...
            FormData(norm_data,datset,class,pos_class,neg_class,i,num_train_fold);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   cross validation   %%%%%%%%%%%%%%%%%%%
        nCenters =3; % not used here
        valid_sz = [10,4,6,8];
        score =zeros(valid_sz);
        margin =zeros(valid_sz);
        error_rate_temp =inf(valid_sz);
        alpha_valid=zeros(valid_sz(1),1);
        NN_valid=zeros(valid_sz(2),1);
        C_valid=zeros(valid_sz(3),1);
        dim_valid =zeros(valid_sz(4),1);
        X=train;
        d_max =min( rank(train), 40);
        d_step = (d_max-1)/8;
        d_step = max(1, d_step);
        d_max = round(d_step*8+1);
        %     dim=SDESVMStruct(k,1).dim;
        for jj = 1:valid_sz(1) % CAUTION:  [jj j c d]=ind2sub([10 4 8 8],ix(1)); 
            
            alpha = jj*0.1;
            alpha_valid(jj) =alpha;
            display(['================  alpha = ' num2str(alpha) ' ================'])
            parfor j=1:valid_sz(2)  %  j=1:4 // [j c d]=ind2sub([5 10 8],ix(1)); % bug? 2012-09-23
                NN=3+j;
                NN_valid(j) =NN;
%                 display(['================  NN = ' num2str(NN) ' ================'])
                [LD] = slle(X', class_train, NN, d_max, alpha);
                sLLE_LD{j,1} = LD;
            end
            for j=1:valid_sz(2)
                NN = NN_valid(j);
                LD = sLLE_LD{j,1};
                training=LD(1:size(train,1),:);
                K_train = training *  training';
                [training, validing] = OutOfSampleExt(K_train, train, valid, group_train, NN,...
                [], d_max, nCenters, [], 'LD', 'LLW', mode);
                for c=1:valid_sz(3)
                    C=2^(2*(c-1));
                    C_valid(c) =C;
%                     display(['================  C = ' num2str(C) ' ================'])
                    parfor d=1:valid_sz(4)
                        dim=round(d_step*d+1);
                        if dim > size(training,2)
                            continue;
                        end
                        dim_valid(d) = dim;
                        LD_train = training(:,1:dim);
                        LD_valid = validing(:,1:dim);
                        [u, gamma, e]=train_svm(LD_train, group_train, C);
                        label_valid =test_svm(LD_train, group_train, LD_valid, u, gamma);
%                         score(jj,j,c,d)=sum(abs(SVMStruct.Alpha))-0.5*SVMStruct.Alpha'*SVMStruct.SupportVectors*SVMStruct.SupportVectors'*SVMStruct.Alpha;
%                         margin(jj,j,c,d)=0.5*SVMStruct.Alpha'*SVMStruct.SupportVectors*SVMStruct.SupportVectors'*SVMStruct.Alpha;
                        error_rate_temp(jj,j,c,d)=sum(group_valid~=label_valid);
                    end
                end
            end
        end
        temp=zeros(size(error_rate_temp));   temp(error_rate_temp~=min(error_rate_temp(:)))=inf;
        ix=find(temp==min(temp(:)));
        [jj j c d]=ind2sub(valid_sz,ix(1)); 
        alpha = alpha_valid(jj); NN=NN_valid(j); C=C_valid(c); dim=dim_valid(d);
        SLLEStruct(k,i).NN=NN;
        SLLEStruct(k,i).dim=dim;
        SLLEStruct(k,i).C=C;
        SLLEStruct(k,i).alpha=alpha;
        display(['================ test ================'])
        save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_sLLE_experiment.mat'],...
            'SLLEStruct','error_rate_sLLE')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%   test               %%%%%%%%%%%%%%%%%%%
        %         display(['================ #' num2str(i) ' ================'])
        %         dim=SDESVMStruct(k,i).dim;
        [LD] = slle(X', class_train, NN, dim, alpha);
        training=LD(1:size(train,1),:);
        K_train = training *  training';
        [training, testing] = OutOfSampleExt(K_train, train, test, group_train, NN,...
            [], dim, nCenters, [], 'LD', method, mode);
        LD_train = training(:,1:dim);
        LD_test = testing(:,1:dim);
        [u, gamma, e]=train_svm(LD_train, group_train, C);
        label_test =test_svm(LD_train, group_train, LD_test, u, gamma);
        error_rate_sLLE(k,i)=sum(group_test~=label_test);
        SLLEStruct(k,i).u =u;
        SLLEStruct(k,i).gamma = gamma;
        SLLEStruct(k,i).MisClassifiedIdx=test_idx(group_test~=label_test,:);
    end
end
sum(error_rate_sLLE,2)
save([name_data '_' num2str(num_samp) 'sample' '_' num2str(num_fold) 'fold_sLLE_experiment.mat'],...
    'SLLEStruct','error_rate_sLLE')
return;
