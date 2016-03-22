% test_WienerSwitchSystem

clear;close all;

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

dbstop if error

opt.nSys = 2; % number of systems
opt.lambda1Init = 0.0001;
opt.lambda1Rate = 10;
opt.epsilon = 0.0;
opt.method = 'moment';
% opt.method = 'convex';
% opt.method = 'convex_noisy';
opt.dataset = 'synthetic';
% opt.dataset = 'mhad';
% opt.dataset = 'msr';

if strcmp(opt.dataset, 'synthetic')
    opt.sysOrders = [2 2]; % order for each system, minimum 2
    opt.numSample = 10; % number of data samples
    opt.switchInd = [5 12 18];
    opt.numDim = 1;
    opt.c1 = 4;
    opt.c2 = 5.09;
    opt.numNeighbors = opt.numSample;
    opt.sysOrd = max(opt.sysOrders);
%     rng('default');
    rng(12);
    [data_clean,r, gt] = switchSysDataGen(opt);
%     data_clean = bsxfun(@minus,data_clean, mean(data_clean,2));
    data_clean_normal = bsxfun(@rdivide,data_clean, max(abs(data_clean),[],2));
    data_linear = data_clean_normal + 0.0 * rand(size(data_clean));
    % data = exp(data_clean)-1 + 0.0 * randn(size(data_clean));
    data = 1./(1+exp(-data_linear)) + opt.epsilon * rand(size(data_clean));
    % data = sqrt(data_clean.^2+5) + 0.0 * randn(size(data_clean));
elseif strcmp(opt.dataset, 'mhad')
    opt.c1 = 0.9;
    opt.c2 = 1.1;
    load ../expData/mhad.mat;
    data = mhad(3,1:opt.numSample);
    data = bsxfun(@minus,data, mean(data,2));
    data = bsxfun(@rdivide,data, max(abs(data),[],2));
    gt = ones(1, opt.numSample);
elseif strcmp(opt.dataset, 'msr')
    opt.sysOrders = [3 3];
    opt.sysOrd = max(opt.sysOrders);
    opt.c1 = 0.9;
    opt.c2 = 1.1;
    [data,gt] = getMSRdata;
end

opt.numNeighbors = size(data, 2);

tic;
if strcmp(opt.method, 'moment')
    opt.lambda1 = 1e-3;
    opt.maxIter = 1000;
    [x,label,rHat,rdHat] = veroneseSDEcvx11(data, opt);
elseif strcmp(opt.method, 'convex')
    opt.lambda1 = 10;
    opt.maxIter = 100;
    [x,label,rHat] = veroneseSDEcvx12(data, opt);
elseif strcmp(opt.method, 'convex_noisy')
    opt.maxIter = 100;
    opt.noiseBound = 0.0001;
    [x,label,rHat] = veroneseSDEcvx13(data, opt);
end
toc

x

figure(2);
plot(label,'x');
hold on;
plot(gt,'o')
hold off;
xlabel('Time step');
ylabel('System Identity');
legend Identified Groundtruth;
title(sprintf('System identification with %s method',opt.method));
