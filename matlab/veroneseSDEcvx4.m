function [x,group] = veroneseSDEcvx4(data, nn, lambda1, lambda2)
% veronese map as constraints

D = pdist2(data',data');
D = D.^2;
% d = size(X, 1);
n = size(data, 2);
Eta = getNNmap(D, nn);
EtaPair = (Eta'*Eta > 0);

epsilon = 1e-6;
d = 3;
p = d*(d+1)/2;
mord = 2;
nSys = 2;
[Dict,Indx] = momentPowers(0,n+p,2*mord);
L = max(Indx);

% build moment matrix
[basis,~] = momentPowers(0,n+p,mord);

Mi = getMomInd(Dict,basis,0,Indx,0);

Ki = getKernelInd4(Dict, n);
% Sh = getStructuredVeroneseMap4(Dict, d, n, 2);
Vi = getVeroneseMap4(Dict, d, n, 2);

h = size(Mi,1);

W1 = eye(h);
W2 = eye(h);
m = zeros(L, 1);
m_pre = ones(L, 1);
terminate = false;
maxIter = 100;
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
%                 ind = getCrossTermIndices3(Dict,i,j);
%                 K(i,j) == m(ind);
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
        sum(K(:)) == 0;
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));
        
%         Vm = reshape(Sh * m, p, []).';
%         sum(Vm, 2) == 0;
        sum(m(Vi)) == 0;

        [P m(Mi); m(Mi)' Q] == semidefinite(2*h);
%         obj =  - trace(K) + lambda1 * norm_nuc(m(Mi)) +  lambda2 * norm_nuc(m(Vi));
        obj = trace(W1*P) + trace(W2*Q);
%         obj = - trace(K) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    s = svd(m(Mi));
    s'
    if s(2)<1e-3*s(1) || norm(m-m_pre)<1e-6 || iter > maxIter
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

% [U,S,V] = svd(K);
% R = S.^0.5 * V';
% s = diag(S);
% c = cumsum(s)/sum(s);
% ind = nnz(c<0.99)+1;
% Y = R(1:ind,:);

x = m(2:n+1);
% Vm = reshape(Sh * m, p, []);
Vm = m(Vi);
X = hankel(m(2:d+1),m(d+1:n+1));
[Vr, powers] = veronese(X,2);
[U,S,V] = svd(Vm);
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