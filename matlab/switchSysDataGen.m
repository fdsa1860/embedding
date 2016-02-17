% switch system data generation

function [data, sys_par] = switchSysDataGen(num_sys, sys_ord, num_sample, num_dim)

%% Data generation

if nargin == 0
    num_sys = 2; % number of systems
    % sys_ord = [2 2 3 3 4 4]; % order for each system, minimum 2
    sys_ord = [2 2]; % order for each system, minimum 2
    num_sample = 5; % number of samples per system
    num_dim = 1;
end
num_Hcol = 10;

sys_par = {};

%% System Generation
theta = (rand(num_sys,1))*2*pi; %
for i = 1:num_sys
    x = rand(sys_ord(i)-2,1);
    p = [cos(theta(i))+1i*sin(theta(i));cos(theta(i))-1i*sin(theta(i));x(:)]; % two complex poles and the rest are real poles
    null{i} = -fliplr(poly([p; rand(num_Hcol-1-sys_ord(i),1)]')); % null space for each system, depends on the number of hankel colomns
    p = -fliplr(poly(p'));
    sys_par{i} = p;
end

%% Generate switch system data
data = zeros(num_sample*num_sys, num_dim);
initValues = rand(sys_ord(1),num_dim) - 0.5;
count = 1;
for i = 1:num_sys
    for j = 1:num_sample
        if i == 1 && j <= sys_ord(1)
            data(count,:) = initValues(j,:);
        else
            data(count,:) = sys_par{i}(1:end-1) * data(count-sys_ord(i):count-1,:);
        end
        count = count + 1;
    end
end
