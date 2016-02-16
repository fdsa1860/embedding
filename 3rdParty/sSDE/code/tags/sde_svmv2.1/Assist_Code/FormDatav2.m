% Partition the dataset into training, testing and validation set.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT:
% norm_data     the normalized data tensor of size D-by-N-by-C, where D is
%               the number of features, N is the number of samples in each
%               class, C is the number of the classes. Note that, when the
%               numbers of samples of classes are not the same, C =1.
% datset        the structure storing the partitions of the dataset.
%               + validation_fold: the index of the validation fold.
%               + fold_idx: the sample index in each fold.
% class         class label of the samples. 
% pos_class     the label for postive class
% neg_class     the label for negative class
% idx           the index of the test fold
% OUTPUT:
% train         the data matrix of the training set N-by-D
% test          the data matrix of the testing set N-by-D
% valid         the data matrix of the validation set N-by-D
% group_train       the group annotation of the training set N-by-1, {+1, -1}
% group_test        the group annotation of the testing set N-by-1, {+1, -1}
% group_valid       the group annotation of the validation set N-by-1, {+1, -1}
% class_train       the class label of the training set N-by-1
% class_test        the class label of the testing set N-by-1
% class_valid       the class label of the validation set N-by-1
% train_idx         the fold index of the training set.
% test_idx          the fold index of the testing set.
% valid_idx         the fold index of the validation set.
function [train, test, valid, group_train, group_test,group_valid, ...
    class_train, class_test, class_valid, train_idx, test_idx, valid_idx]= ...
    FormDatav2(norm_data,datset,class,pos_class,neg_class,idx)
train=[];test=[];valid=[];
group_train=[];group_test=[];group_valid=[];
class_train=[];class_test=[];class_valid=[];
train_idx=[];test_idx=[];valid_idx=[];

if length(size(norm_data))==2
    pos_class=1; neg_class=[];
end
valid_idx=datset(pos_class).validation_fold;
test_idx=datset(pos_class).fold_idx(:,idx);
train_idx=datset(pos_class).fold_idx;
train_idx(:,idx)=[];
train_idx=train_idx(:);    test_idx= test_idx(:);    valid_idx=valid_idx(:);

train=squeeze(norm_data(:,train_idx,pos_class)'); test=squeeze(norm_data(:,test_idx,pos_class)');
valid=squeeze(norm_data(:,valid_idx,pos_class)');
class_train=class(train_idx,pos_class);         class_test=class(test_idx,pos_class);
class_valid=class(valid_idx,pos_class);

train_idx=[ones(size( train_idx(:)))*pos_class,train_idx(:)]; 
test_idx=[ones(size( test_idx(:)))*pos_class,test_idx(:)]; 
valid_idx=[ones(size( valid_idx(:)))*pos_class,valid_idx(:)]; 

if ~isempty(neg_class)
    valid_neg_idx=datset(neg_class).validation_fold;
    test_neg_idx=datset(neg_class).fold_idx(:,idx);
    train_neg_idx=datset(neg_class).fold_idx;
    train_neg_idx(:,idx)=[];
    
    train_neg_idx=train_neg_idx(:);    test_neg_idx= test_neg_idx(:);    valid_neg_idx=valid_neg_idx(:);
    train=[train; squeeze(norm_data(:,train_neg_idx,neg_class))'];
    test=[test; squeeze(norm_data(:,test_neg_idx,neg_class))'];
    valid=[valid; squeeze(norm_data(:,valid_neg_idx,neg_class))'];
    class_train=[class_train; class(train_neg_idx,neg_class)];
    class_test=[class_test; class(test_neg_idx,neg_class)];
    class_valid=[class_valid; class(valid_neg_idx,neg_class)];
    
    train_idx=[train_idx;[ones(size(train_neg_idx))*neg_class, train_neg_idx]];
    test_idx=[test_idx;[ones(size(test_neg_idx))*neg_class, test_neg_idx]];
    valid_idx=[valid_idx;[ones(size( valid_neg_idx))*neg_class, valid_neg_idx]];
end
group_train=-ones(size(class_train));    group_train(class_train==pos_class)=1;
group_test=-ones(size(class_test));    group_test(class_test==pos_class)=1;
group_valid=-ones(size(class_valid));    group_valid(class_valid==pos_class)=1;
return;