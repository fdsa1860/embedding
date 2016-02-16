% test_WienerSwitchSystem

clear;close all;

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

dbstop if error

nSys = 2; % number of systems
% sys_ord = [2 2 3 3 4 4]; % order for each system, minimum 2
sysOrd = [2 2]; % order for each system, minimum 2
numSample = 5; % number of samples per system

% rng('default');
rng(1);
[data_clean,r] = switchSysDataGen(nSys, sysOrd, numSample);
numNeighbors = 8;

data = 10*data_clean + 0.0 * randn(size(data_clean));
% data = exp(data_clean)-1 + 0.0 * randn(size(data_clean));
% data = (data_clean*10).^3 + 0.0 * randn(size(data_clean));

% data = data - mean(data);
figure(1);
plot(data_clean,'*');
hold on;
plot(data,'g*');
hold off;

% maxSysOrd = max(sysOrd);
% H = hankel_mo(data', [maxSysOrd+1 length(data)-maxSysOrd]);
% D = pdist2(H', H');
D = pdist2(data,data);

tic;
lambda1 = 100;
lambda2 = 2;
% [Y,group] = veroneseSDEcvx(D, numNeighbors, lambda1, lambda2);
[x,group] = veroneseSDEcvx5(D, numNeighbors, lambda1, lambda2);
toc

x
figure(1);hold on;
plot(x(1,:),'-*');
hold off;

figure(2);
plot(group,'x');