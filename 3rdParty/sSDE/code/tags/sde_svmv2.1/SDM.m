% SVDM function. Both training and testing on one set of train, test and
% validation set partition.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
% train     the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% test      the high dimension test data, which is Ng-by-D, where Ng is the number
%           of samples
% class     the class annotation of the samples.
% dim: the dimensionality of the projected space
% D: the weight of the hinge loss term 
% mu: the target classification margin.
%   
% OUTPUT
% theta: sign(Z_train*theta) is the predicted label 
% Z_train, W: W*Z_train is used to approximated the data matrix 
% label_train: NOT USED
% label_test: the predicted label annotation of the test set.
function [W, Z_train, theta, label_train,label_test]=SDM(train, test,class, dim, D, mu)
label_train=[];label_test=[];
[groupIndex, groupString] = grp2idx(class);
num_class=length(groupString);
[W,Z_train,theta]=SDM_train(train, class, dim, num_class, D, mu);
Z_test=(W'\test')';
temp=Z_test*theta;
temp(temp<0)=0;
[temp ix]=sort(temp,2,'descend');
label_test=zeros(size(Z_test,1),1);
label_test(ix(:,1)==1)=1;
label_test(ix(:,1)==2)=-1;
return;

% SVDM training function.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT
% X         the high dimension trainning data, which is Nt-by-D, where Nt is the number
%           of samples
% class     the class annotation of the samples
% num_class the number of class
% dim: the dimensionality of the projected space
% D: the weight of the hinge loss term 
% mu: the target classification margin.  
% OUTPUT
% theta_old: sign(Z_train*theta) is the predicted label 
% Z_old, W_old: W*Z_train is used to approximated the data matrix 
function [W_old,Z_old,theta_old]=SDM_train(X, class, dim, num_class, D, mu)
[num_smp,num_dim]=size(X);
%% initialization
Y=-1*ones(num_smp,num_class);
% for i=1:num_class
%     Y(class==i,i)=1;
% end
Y(class==1,1)=1;
Y(class==-1,2)=1;
temp=randn(dim+1,num_class);
temp2=diag(1./sqrt(diag(temp'*temp)));
theta=temp*temp2;

W=randn(dim+1,num_dim);

Z=ones(num_smp,dim+1);
temp=randn(num_smp,dim);
temp2=diag(1./sqrt(diag(temp*temp')));
temp=temp2*temp;
Z(:,2:dim+1)=temp;

%%
Z_old=Z+1;
k=0;
while (k<2 || val(k-1)-val(k)>1e-1)
%     display(num2str(val));
%     norm(Z_old-Z,'fro')
    %% solve W
    W_old=W;
    for i=1:num_dim
        W(:,i)=Z\X(:,i);
    end
    
    %% solve theta
    theta_old=theta;
    % theta=zeros(num_smp,num_class);
    for j=1:num_class
        cvx_begin sdp quiet
        cvx_solver  sedumi %sdpt3 %sedumi
        cvx_precision medium
        variable h(num_smp,1);
        variable t(dim+1,1)
        expression rho(num_smp,1);
        rho=diag(Y(:,j))*Z*t;
        
        minimize (sum(h))
        
        subject to
        h>=0;
        h>=D*(mu-rho);
        t'*t<=1;
        
        cvx_end
        theta(:,j)=t;
    end
    %% solve Z
    Z_old=Z;
    for i=1:num_smp
%         i
        cvx_begin sdp quiet
        cvx_solver sedumi
        variable z(1,dim+1);
        variable h(1,num_class);
        expression rho(1,num_class);
        rho=Y(i,:).*(z*theta(:,j));
        minimize ((X(i,:)-z*W)*(X(i,:)-z*W)'+sum(h))
        subject to
        h>=0;
        h>=D*(mu-rho);
        z(1)==1;
        z(2:dim)*z(2:dim)'<=1;
        cvx_end
        Z(i,:)=z;
    end
    % % for i=1:num_smp
    %     cvx_begin
    %     cvx_solver sedumi
    %     variable z(num_smp,dim);
    %     variable h(num_smp,num_class);
    %
    %     expression rho(num_smp,num_class);
    %     rho=Y.*([ones(num_smp,1) z]*theta);
    % %     minimize (trace((X-[ones(num_smp,1) z]*W)*(X-[ones(num_smp,1) z]*W)')+sum(h(:)))
    %     minimize (norm((X-[ones(num_smp,1) z]*W),'fro')+sum(h(:)))
    %     subject to
    %     h>=0;
    %     h>=D*(mu-rho);
    %     diag(z*z')<=1;
    %     cvx_end
    % % end
    k=k+1;
    rho=Y.*(Z*theta);
    h=D*(mu-rho);
    h(h<0)=0;
    val(k)=norm(X-Z*W,'fro')^2+sum(h(:));
end
return;