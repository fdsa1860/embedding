% Form the local isometric constraints 
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUTs:
%       eta: the neighborhood matrix. eta_ij =1 if and only if there is an
%           edge between i-th and j-th sample.
%   
% OUTPUTs
%       S_eta: coefficient vectors for the local isometric constraints. So
%           that for each edge (eta_ij =1) there is a row in S_eta that the
%           i-th and j-th entry in that row are 1 and -1 respectively and
%           all the other entries are 0.
%       relax: the relaxation coefficients of the local geometric
%               constraints. NOT USED
function [S_eta, relax]=FormConstraintAll(eta,varargin)
    if size(varargin)>0
        Y=varargin{1};
        ratio=varargin{2};
    end
    eta=eta|eta';
    eta=eta-diag(diag(eta));
    eta=triu(eta);
    
    N_eta=sum(eta(:));
    [idxx idxy]=ind2sub(size(eta),find(eta(:)));
    temp=sparse(1:length(idxx),idxx,ones(size(idxx)),N_eta, size(eta,1));
    S_eta=temp+sparse(1:length(idxy),idxy,-ones(size(idxy)),N_eta, size(eta,1));
    
    relax=[];
%     temp=S_eta*Y;
%     relax=ones(size(temp));
%     relax(temp~=0,:)=1+ratio;
return;