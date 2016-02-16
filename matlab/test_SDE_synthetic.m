% try SDE on synthetic data of different systems

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

rng('default');
data_generation;
data_clean = data;
n = 50;
d = num_frame;
N = n * num_sys;
nc = 30;
% opt.metric = 'binlong';
% opt.metric = 'AIRM';
% opt.metric = 'LERM';
opt.metric = 'JBLD';
% opt.metric = 'KLDM';
opt.sigma = 1e-2;

numNeighbors = 8;

e = 1;
data = data_clean + e * noise_data;

HH = cell(1,size(data,1));
for i = 1:size(data,1)
    H1 = hankel_mo(data(i,:), [d-nc+1, nc]);
    H1_p = H1 / (norm(H1*H1','fro')^0.5);
    HH1 = H1_p' * H1_p;
    if strcmp(opt.metric,'AIRM') || strcmp(opt.metric,'LERM')...
            || strcmp(opt.metric,'KLDM') || strcmp(opt.metric,'JBLD')
        HH{i} = HH1 + opt.sigma(1) * eye(nc);
    elseif strcmp(opt.metric,'binlong')
        HH{i} = HH1;
    end
end

D = HHdist(HH,[],opt);

tic;
Y = mySDEcvx(D, numNeighbors);
toc

% SVM classification
C_val = 10;
total_accuracy = zeros(1, num_fold);
cw_accuracy = zeros(num_fold, 6);
confusion_matrices = cell(1, num_fold);
for iFold = 1:num_fold
    
    trainData = [];
    testData = [];
    y_train = [];
    y_test = [];
    for i = 1:num_sys
        trainData = [trainData Y(:,index_train(:,iFold)+(i-1)*n)];
        y_train = [y_train label(index_train(:,iFold)+(i-1)*n)'];
        testData = [testData Y(:,index_test(:,iFold)+(i-1)*n)];
        y_test = [y_test label(index_test(:,iFold)+(i-1)*n)'];
    end
    
    [total_accuracy(iFold), cw_accuracy(iFold,:), confusion_matrices{iFold}] =...
    svm_one_vs_all(trainData, testData, y_train', y_test', C_val);
end
disp('average total accuracy is ');
mean(total_accuracy)

