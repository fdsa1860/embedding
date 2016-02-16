% Function/Wrapper to train with matlab SVM command.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUT:
% Training  the training data representatives in the output space, each
%           column is a embedded data sample.
% y         the label of each sample {-1 +1}.
% C         the weight of the SVM error term
% OUTPUT:
% u         the coefficient for support vectors.
% gamma     the offset of the boundary.
% e         the errors.
% score     the objective value.
% margin    the objective on margin.
function [u, gamma, e, score, margin]=train_svm(Training, y,C)
D=diag(y);
class=y; class(class==-1)=2;
u=zeros(size(y));
% train SVM
options = statset('MaxIter',1e7);
SVMStruct = svmtrain(Training,class,'kernel_function','linear',...
    'boxconstraint',C, 'autoscale', false,'kktviolationlevel', 0.05,...
    'tolkkt', 5e-3, 'options',options);
u(SVMStruct.SupportVectorIndices,1)=SVMStruct.Alpha./y(SVMStruct.SupportVectorIndices);
gamma=-SVMStruct.Bias;
e=1-D*(Training*Training'*D*u-gamma);
e(e<=0)=0;

score=sum(abs(SVMStruct.Alpha))-0.5*SVMStruct.Alpha'*SVMStruct.SupportVectors*SVMStruct.SupportVectors'*SVMStruct.Alpha;
margin=0.5*SVMStruct.Alpha'*SVMStruct.SupportVectors*SVMStruct.SupportVectors'*SVMStruct.Alpha;

% [u gamma e label] = KSVMv1 (K(1:Nt,1:Nt),y, lam2);
return;

% %% train Lib-SVM
% function [u, gamma, e, score, margin]=train_svm(Training, y,C)
% e=[];
% D=diag(y);
% class=y; class(class==-1)=2;
% u=zeros(size(y));
% % train SVM
% options=['-s 0 -t 0 -e 0.0001 -c ' num2str(C) ' -m 1000'];
% model =libsvmtrain([], y, Training, options);
% if model.Label(1)>0
%     w=model.sv_coef'*model.SVs;
%     gamma= -model.rho;
% else
%     w=-model.sv_coef'*model.SVs;
%     gamma= model.rho;
% end
% w=w';
% % [u gamma e label] = KSVMv1 (K(1:Nt,1:Nt),y, lam2);
% return;