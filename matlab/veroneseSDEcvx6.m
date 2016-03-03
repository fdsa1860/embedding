function [x,group] = veroneseSDEcvx6(data, nn, lambda1, lambda2)
% solve 2 systems switch system with veronese map and moment

D = pdist2(data',data');
D = D.^2;
% d = size(X, 1);
n = size(data, 2);
Eta = getNNmap(D, nn);
EtaPair = (Eta'*Eta > 0);

epsilon = 1e-6;
d = 3;

q = n*(n+1)/2; % number of kernel variables
p = d*(d+1)/2; % number of veronese map variables

% build moment matrix
mord = 1;
nSys = 2;
[Dict,Indx] = momentPowers(0,q+p,2*mord);
L = max(Indx);
[basis,~] = momentPowers(0,q+p,mord);
Mi = getMomInd(Dict,basis,0,Indx,0);
h = size(Mi,1);

% get indices for regressors
rInd = [1 q+2:q+p+1];
Ri = Mi(rInd, rInd);
% h = size(Ri,1);

% get Kernel indices
Ki = getKernelInd(n);

Vi = getStructuredVeroneseMap6(Dict, d, n, 2);



W1 = eye(h);
W2 = eye(h);
m = zeros(L, 1);
m_pre = ones(L, 1);
terminate = false;
maxIter = 10;
iter = 1;
while ~terminate
    m_pre = m;
    cvx_begin sdp
    cvx_solver mosek
        variables m(L, 1);
        variables P(h, h) Q(h, h);
        K = m(Ki);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
%                     K(i,i)+K(j,j)-2*K(i,j) == D(i,j);
%                     K(i,i)+K(j,j)-2*K(i,j) >= D(i,j)-epsilon;
%                     K(i,i)+K(j,j)-2*K(i,j) <= D(i,j)+epsilon;
                    K(i,i)+K(j,j)-2*K(i,j) <= 25*D(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) >= 16*D(i,j);
                end
            end
        end
        K == semidefinite(n);
%         sum(K(:)) == 0;
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));
        
        Vm = m(Vi);
        sum(Vm, 2) == 0;
%         sum(Vm, 2) >= -1e-6;
%         sum(m(q+2:q+p+1)) == 1;
        m(q+2) == 1;
%         ind = find(sum(Dict,2)==2 & sum(Dict(:,1:q)==2,2));
%         m(ind) <= 1;

        [P m(Mi); m(Mi)' Q] == semidefinite(2*h);
        
%         obj =  - trace(K) + lambda1 * norm_nuc(m(Mi)) +  lambda2 * norm_nuc(m(Vi));
        obj = trace(W1*P) + trace(W2*Q) ;
%         obj = sum(m(ind)) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    s = svd(m(Mi));
    s'
    if s(2)<1e-5*s(1) || norm(m-m_pre,inf)<1e-5 || iter > maxIter
        terminate = true;
    end
    iter = iter + 1;
    
    sy = svd(P);
    sz = svd(Q);
    W1 = inv(P + sy(2)*eye(h));
    W2 = inv(Q + sz(2)*eye(h));
    scale = norm(blkdiag(W1,W2),2);
    W1 = W1/scale;
    W2 = W2/scale;
end

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
X = R(1:ind,:);

x = R(1,:);
X = hankel(x(1:d),x(d:end));
[Vr, powers] = veronese(X,2);
[U,S,V] = svd(Vr);
s = diag(S)
% assert(s(end)<1e-6);
r = U(:,end);

%compute only one normal at each point to the subspace containing the point,
%assuming the codim is one.
% Dpn=derivative(c,powers,x);
[Dpn,normDpn] = cnormalize(derivative(r,powers,X));

method = 'Cos^2';
%compute the similarity matrix
switch method
    case 'Cos'
        affMat=(abs(Dpn'*Dpn));
    case 'Cos^2'
        affMat=(abs(Dpn'*Dpn)).^2;
    case 'Exp_-sin^2'
        affMat=exp((abs(Dpn'*Dpn)).^2-1); % i.e., exp(-sin^2)
end

%segment using spectral clustering
group=SpectralClustering(affMat,2);

end