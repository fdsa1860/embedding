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
    FormData(norm_data,datset,class,pos_class,neg_class,idx,varargin)
num_train_fold=[];
if size(varargin)>0
    num_train_fold=varargin{1};
end
if  ~isempty(num_train_fold) && num_train_fold>1
idx=[idx:idx+num_train_fold-1]-1;
idx=mod(idx,datset.num_fold)+1;
end
[train, test, valid, group_train, group_test,group_valid, ...
    class_train, class_test, class_valid, train_idx, test_idx, valid_idx]= ...
    FormDatav2(norm_data,datset,class,pos_class,neg_class,idx);
