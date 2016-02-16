% Function of testing samples with trained SVM.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT:
% Training  the training data representatives in the output space, each
%           column is a embedded data sample.
% y         the label of each sample {-1 +1}.
% Testing   the testing data representatives in the output space, each
%           column is a embedded data sample.
% u         the coefficient for support vectors.
% gamma     the offset of the boundary.
% OUTPUT:
% label_test    the predicted label for the testing data, {+1 -1}.
function [ label_test ]=test_svm(Training, y, Testing, u, gamma)
label_test = sign(Testing *Training' *diag(y) *u - gamma);
return;
% %% test Lib-SVM
% function [ label_test ]=test_svm(Training, y, Testing, w, gamma)
% label_test = sign(Testing*w + gamma);
% return;