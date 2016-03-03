% test_WienerSwitchSystem

clear;close all;

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

dbstop if error

nSys = 2; % number of systems
% sys_ord = [2 2 3 3 4 4]; % order for each system, minimum 2
sysOrd = [2 2]; % order for each system, minimum 2
numSample = 5; % number of samples per system
numDim = 3;

% rng('default');
rng(1);
[data_clean,r] = switchSysDataGen(nSys, sysOrd, numSample, numDim);
% data_clean = 10 * data_clean;
data_clean = bsxfun(@minus,data_clean, mean(data_clean));
data_clean = bsxfun(@rdivide,data_clean, max(abs(data_clean)));
data_clean = data_clean.';
numNeighbors = 9;

% data = data_clean + 0.0 * rand(size(data_clean));
% data = exp(data_clean)-1 + 0.0 * randn(size(data_clean));
data = 1./(1+exp(-data_clean)) + 0.0 * randn(size(data_clean));
% data = sqrt(data_clean.^2+5) + 0.0 * randn(size(data_clean));

% data = data - mean(data);
figure(1);
plot(data_clean(1,:),'*');
hold on;
plot(data(1,:),'g*');
hold off;

tic;
lambda1 = 1000;
lambda2 = 2;
% [Y,group] = veroneseSDEcvx(D, numNeighbors, lambda1, lambda2);
[x,group] = veroneseSDEcvx11(data, numNeighbors, lambda1, lambda2);
toc

x
figure(1);hold on;
plot(x(1,:),'-*');
hold off;

figure(2);
plot(group,'x');