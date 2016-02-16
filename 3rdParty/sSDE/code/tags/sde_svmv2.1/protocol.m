% script to set up the experiment protocol

% function [exp_pair validation_fold test_fold]= protocol(name_data,num_fold)
switch name_data
    case 'usps'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   usps               %%%%%%%%%%%%%%%%%%%
% cross validation protocle
        exp_pair=[1 2; 1 3; 2 8; 8 9];
        validation_fold=1;
        test_fold=2:num_fold;
    case 'PAL'
        exp_pair=[1 2];
        validation_fold=1;
        test_fold=2:num_fold;
    case 'Ionosphere'
        exp_pair=[1 2];
        validation_fold=1;
        test_fold=2:num_fold;
    otherwise
        display(name_data);
        exp_pair=[1 2];
        validation_fold=1;
        test_fold=2:num_fold;
end
% return;