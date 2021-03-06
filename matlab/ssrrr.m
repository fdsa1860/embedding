function [label, R1, R2] = ssrrr(K,ord,epsilon)

% load('K');
% K = K(1:20,1:20);
m = size(K,2);

if nargin == 1
    ord = 2;    % define system order
    epsilon = 0.6; % set noise bound
end

% form data
N = m-ord;
X = cell(1,N);
for i = 1:N
    X{i} = K(:,i:i+ord);
end

label = zeros(N,1);

% run first greedy iteration
disp('Begin to run the first greedy iteration.');
w = ones(1,N);
rpre = zeros(ord+1,N);
for iter = 1:20
    cvx_begin
    cvx_solver gurobi;
    variable R(ord+1,1);
    variable r(ord+1,N);
    expression F(N,1);
    
    for i = 1:N
        %         F(i) = max(abs(r(:,i)-R));
        F(i) = norm(r(:,i)-R);
    end
    
    minimize w*F
    subject to
    R(1) == 1;
    for i = 1:N
        abs(X{i}*r(:,i))<=epsilon;
    end
    cvx_end
    
    if norm(rpre(:)-r(:),inf)<=1e-4
        break;
    end
    rpre = r;
    
    w = 1./(F+1e-2)';
end

% record regressor in first greedy iteration
R1 = R;

% remove the first greedy cluster
del = F<=1e-3;
label(del) = 1;

X(del) = [];
N = numel(X);

% run second greedy iteration
disp('Begin to run the second greedy iteration.');
w = ones(1,N);
rpre = zeros(ord+1,N);
for iter = 1:20
    cvx_begin
    cvx_solver gurobi;
    variable R(ord+1,1);
    variable r(ord+1,N);
    expression F(N,1);
    
    for i = 1:N
        %         F(i) = max(abs(r(:,i)-R));
        F(i) = norm(r(:,i)-R);
    end
    
    minimize w*F
    subject to
    R(1) == 1;
    for i = 1:N
        abs(X{i}*r(:,i))<=epsilon;
    end
    cvx_end
    
    if norm(rpre(:)-r(:),inf)<=1e-4
        break;
    end
    rpre = r;
    
    w = 1./(F+1e-2)';
end

% record regressor in second greedy iteration
R2 = R;

% remove the second greedy cluster
del = F<=1e-3;
t = zeros(1,N);
t(del) = 2;
label(label==0) = t;

X(del) = [];

end