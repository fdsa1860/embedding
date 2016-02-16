clear all; clc;
colordef none;
colormap jet;
figure(1); clf; drawnow;
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off');

% LLE 
N=1000; 
D=3;
d=2;
K=8; 
tol=0.01;
markerSize=12;
rand('state',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S MANIFOLD

X = zeros(D,N);
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
angle = pi*(1.5*rand(1,N/2)-1);
height = 5*rand(1,N);
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

% LLE
fprintf('S-MANIFOLD: K=%d, N=%d, tol=%g\n',K,N,tol);
[Y,eigenvals,neighbors] = lle(X,K,d,tol);

% SCATTERPLOT
subplot(3,3,3); cla; hold on;
scatter(Y(1,:),Y(2,:),markerSize,color,'o','filled');
view(2); grid off; axis tight; drawnow;
title('(C)','FontSize',18);
set(gca,'FontSize',18);
set(gcf,'PaperPositionMode','auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWIN PEAKS

inc = 0.1;
[xx2,yy2] = meshgrid(-1:inc:1);
zz2 = sin(pi*xx2).*tanh(3*yy2);
xy = 1-2*rand(2,N);
xx = [xy; sin(pi*xy(1,:)).*tanh(3*xy(2,:))];
color = xx(3,:);
subplot(3,3,4); cla; 
surf(xx2,yy2,zz2); 
grid off; 
view(3); 
set(gca,'FontSize',18);
drawnow;

subplot(3,3,5); cla; hold on;
scatter3(xx(1,:),xx(2,:),xx(3,:),markerSize,color,'o','filled'); 
% scatter3(xx(1,:),xx(2,:),xx(3,:),markerSize,'ko');
view(3); 
grid off;
set(gca,'FontSize',18);
drawnow;

% LLE
fprintf('TWIN PEAKS: K=%d, N=%d, tol=%g\n',K,N,tol);
[Y,eigenvals,neighbors] = lle(xx,K,d,tol);

% SCATTERPLOT
subplot(3,3,6); cla; hold on;
scatter(Y(1,:),Y(2,:),markerSize,color,'o','filled');
view(2); grid off; axis tight; drawnow;
set(gca,'FontSize',18);
set(gcf,'PaperPositionMode','auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUNCTURED SPHERE

inc = 1/4;
[xx,yy] = meshgrid(-5:inc:5);
rr2 = xx(:).^2 + yy(:).^2;
[tmp ii] = sort(rr2);
Y = [xx(ii(1:N))'; yy(ii(1:N))'];

subplot(3,3,7); cla;
[xx,yy,zz] = sphere(20);
zz=zz+1;
surf(xx(1:16,:),yy(1:16,:),zz(1:16,:));
axis equal;
grid off;
set(gca,'ZLim',[0 2]);
set(gca,'FontSize',18);
drawnow;

% STEREOGRAPHIC MAPPING
subplot(3,3,8); cla; 
a = 4./(4+sum(Y.^2));
X = [a.*Y(1,:); a.*Y(2,:); 2*(1-a)];
color = X(3,:);
scatter3(X(1,:),X(2,:),X(3,:),markerSize,color,'o','filled'); hold on;
axis equal;
grid off;
set(gca,'FontSize',18);
set(gca,'ZLim',[0 2]);
drawnow;

% LLE
fprintf('PUNCTURED SPHERE: K=%d, N=%d, tol=%g\n',K,N,tol);
[Y,eigenvals,neighbors] = lle(X,K,d,tol);
subplot(3,3,9); cla;
scatter(Y(1,:),Y(2,:),markerSize,color,'o','filled'); hold on;
set(gca,'FontSize',18);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p=1:9
  subplot(3,3,p);
  set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
end;
set(gcf,'Position',[320 90 626 628]);
% print -djpeg manifolds.jpg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
