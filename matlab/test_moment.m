% test_moment

clear;close all;

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

rng('default');

mord = 2;
[Dict,Indx] = momentPowers(0,2,2*mord);
L = max(Indx);

% build moment matrix
[basis,~] = momentPowers(0,2,mord);

Mi = getMomInd(Dict,basis,0,Indx,0);

% build equality constraint x^2=4 (x^2-4>=0 and x^2-4<=0)
[local,~] = momentPowers(0,2,mord-1);
powerlist_local = [0,0;2,0];
coef_local = [-4,1];
Li = cell(2,1);
for i = 1:size(powerlist_local,1)
    Li{i} = getMomInd(Dict,local,powerlist_local(i,:),Indx,0);
end
    

% build the objective function
powerlist_obj = [2,0;1,1;0,2];
coef_obj = [1,2,1];
obj = zeros(1,L);

for i = 1:size(powerlist_obj,1)
    j = getVecInd(Dict,powerlist_obj(i,:),Indx,0);
    obj(j) = coef_obj(i);
end

% solve the moment SDP
cvx_begin sdp
cvx_solver sedumi;
variable m(L,1);

minimize obj*m
subject to
m(Mi(1,1))==1;
m(Mi)>=0;
coef_local(1)*m(Li{1})+coef_local(2)*m(Li{2})==0;
cvx_end

