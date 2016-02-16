% toy example of SDE
% Input:
% X: d by n matrix
% Output:
% Y: r by n matrix

addpath(genpath('../3rdParty'));

colordef none;
colormap jet;
figure(1); clf; drawnow;
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off');

% LLE 
n=200; 
D=3;
d=2;
K=8; 
tol=0.01;
markerSize=12;
rand('state',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S MANIFOLD

X = zeros(D,n);
colordef none; 
colormap jet;

% MANIFOLD
tt = [-1:0.1:0.5]*pi;
uu = tt(end:-1:1);
hh = [0:0.1:1]*5;
xx = [cos(tt) -cos(uu)]'*ones(size(hh));
zz = [sin(tt) 2-sin(uu)]'*ones(size(hh));
yy = ones(size([tt uu]))'*hh;
cc = [tt uu]' * ones(size(hh));
subplot(3,3,1); cla; 
surf(xx,yy,zz,cc);
view([12 -20 3]); grid off; 
set(gca,'ZTick',[-1:3]);
set(gca,'FontSize',18);
title('(A)','FontSize',18);
drawnow;

% DATA
angle = pi*(1.5*rand(1,n/2)-1);
height = 5*rand(1,n);
color = [angle angle];
X(1,:) = [cos(angle) -cos(angle)];
X(2,:) = height; 
X(3,:) = [sin(angle) 2-sin(angle)];

% SCATTERPLOT
subplot(3,3,2); cla;
scatter3(X(1,:),X(2,:),X(3,:),markerSize,color,'o','filled'); hold on;
view([12 -20 3]); grid off; 
set(gca,'ZTick',[-1:3]);
set(gca,'FontSize',18);
title('(B)','FontSize',18);
drawnow;

% % LLE
% fprintf('S-MANIFOLD: K=%d, N=%d, tol=%g\n',K,N,tol);
% [Y,eigenvals,neighbors] = lle(X,K,d,tol);

% tic;
% [LD, optvalcvx] = SDEcvx(X', 1e-6, 8);
% Y = LD';
% toc

tic;
D = pdist2(X',X');
Y = mySDEcvx(D,8);
toc

% % d = 8;
% % n = 100;
% % X = rand(d, n);
% G = X' * X;
% g = diag(G);
% G_ii = kron(ones(n,1), g(1:n));
% G_jj = kron(g(1:n), ones(n,1));
% G_ij = G(:);
% 
% cvx_solver sedumi
% cvx_begin
%     variables K(n, n);
%     K == semidefinite(n);
%     obj = - trace(K);
%     sum(K(:)) == 0;
%     k = diag(K);
%     K_ii = kron(ones(n,1), k(1:n));
%     K_jj = kron(k(1:n), ones(n,1));
%     K_ij = K(:);
%     K_ii + K_jj - 2*K_ij == G_ii + G_jj - 2*G_ij;
%     minimize(obj);
% cvx_end
% [U,S,V] = svd(K);
% R = S.^0.5 * V';
% s = diag(S);
% c = cumsum(s)/sum(s);
% ind = nnz(c<0.99)+1;
% Y = R(1:ind,:);

% SCATTERPLOT
subplot(3,3,3); cla; hold on;
scatter(Y(1,:),Y(2,:),markerSize,color,'o','filled');
view(2); grid off; axis tight; drawnow;
title('(C)','FontSize',18);
set(gca,'FontSize',18);
set(gcf,'PaperPositionMode','auto');
