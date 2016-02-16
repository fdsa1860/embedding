% try SDE on synthetic data of different systems

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

numNeighbors = 10;

rng('default');
r = 1;
theta = 0:0.1:2*pi;
x = cos(theta);
y = sin(theta);
plot(x,y,'*');
z = zeros(1, length(theta));
s = rand(1, length(theta));
z(s>0.5) = 1;
% z(s<=0.5) = -1;
figure(1);
plot3(x,y,z,'-*');

data_clean = [x; y; z];
data = data_clean + 0.0 * randn(size(data_clean));
D = pdist2(data', data');

tic;
Y = lowRankSDEcvx(D, numNeighbors, 1);
toc

figure(2);
plot3(Y(1,:),Y(2,:),zeros(1, length(theta)),'-*');