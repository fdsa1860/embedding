% Semi-Definite Embedding with CVX toolbox
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% INPUTs:
% S1        the high dimension data, which is N-by-D, where N is the number
%           of samples
% NN        the number of local neightbors.
% eps       relax the local isometric constraint
%
% OUTPUTs
% LD        Low dimension mapped data.

function [LD optvalcvx] = SDEcvx(S1, eps, NN)
% initialization
S1=double(S1);
[N D]=size(S1);
LD=[];
%% Rank Minimization with Kernel Matrix
% S1=S1-repmat(mean(S1,1),N,1);
G=S1*S1';
% normalize the S1's energy (2-norm) to N.
% S1=S1./sqrt(trace(G)/N);
% G=G/trace(G)*N;

% Neighborhood Matrix
W=NNeighbour(S1',NN);
eta=W';
eta=double(abs(eta)>0);
[eta]=FullConn(S1, eta);
eta=(eta>0)|((eta'*eta)>0);
eta=eta|eta';
%%%%%%%%%%
[MC ,~]=FormConstraintAll(eta);
MGM=sum((MC*G).*MC,2);

cvx_begin sdp quiet
cvx_solver   sedumi
cvx_precision medium
variable K(N,N) symmetric;
minimize (-trace(K))

subject to
K==semidefinite(N);
sum(sum(K))==0;
sum((MC*K).*MC,2)<=(1+eps)*MGM;
if eps <1
    sum((MC*K).*MC,2)>=(1-eps)*MGM;
end
cvx_end
disp(cvx_status)

if (strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) %&& trace(K)>N
    optvalcvx=cvx_optval;
    [U D V]=svd(K);
    LD=U*D.^0.5;
else
    return
end


return;