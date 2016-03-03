function [x,group] = veroneseSDEcvx5(data, k, lambda1, lambda2)
% do not use moment, detect switching between two systems by minimizing
% rank of veronese matrix

D = pdist2(data',data');
D2 = D.^2;
% d = size(X, 1);
n = size(data, 2);
Eta = getNNmap(D, k);
EtaPair = (Eta'*Eta > 0);

epsilon = 1e-6;
d = 3;
p = d*(d+1)/2;

% Sh = getStructuredVeroneseMap5(n, d, 2);
Vi = getVeroneseMap5(n, d);

W1 = eye(n-d+1);
W2 = eye(p);
WK = eye(n);
K = zeros(n, n);
K_pre = zeros(n,n);
maxIter = 10;
iter = 1;
terminate = false;
while ~terminate
    cvx_begin sdp
    cvx_solver mosek
        variables K(n, n);
        variables P(n-d+1, n-d+1) Q(p, p);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
%                     K(i,i)+K(j,j)-2*K(i,j) == D2(i,j);
%                     K(i,i)+K(j,j)-2*K(i,j) >= D2(i,j)-epsilon;
%                     K(i,i)+K(j,j)-2*K(i,j) <= D2(i,j)+epsilon;
                    K(i,i)+K(j,j)-2*K(i,j) <= 25*D2(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) >= 16*D2(i,j);
                end
            end
        end
        K == semidefinite(n);
%         sum(K(:)) == 0;
%         Vr = reshape(Sh * K(:), p, n-d+1).';
        Vr = K(Vi);

        [P Vr; Vr' Q] == semidefinite(p+n-d+1);
        
%         norm_K3 = 0;
%         for i=1:n-d+1
%             norm_K3 = norm_K3 + norm_nuc(K(i:i+d-1,i:i+d-1));
%         end
%         obj =  norm_K3 + lambda1 * (trace(W1*P) + trace(W2*Q));
%         obj = (trace(W1*P) + trace(W2*Q)) ;
%         obj = - trace(K) + lambda1 * (trace(W1*P) + trace(W2*Q));
        obj = trace(WK*K) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    s = svd(Vr);
    s'
    sk = svd(K);
    sk'

    
    sp = svd(P);
    sq = svd(Q);
    W1 = inv(P + sp(5)*eye(n-d+1));
    W2 = inv(Q + sq(5)*eye(p));
    scale = norm(blkdiag(W1,W2),2);
    W1 = W1/scale;
    W2 = W2/scale;
    
    WK = inv(K + sk(2)*eye(n));
    WK = WK / norm(WK);
    
    iter = iter + 1;
    
    %     if s(end)/sum(s) < 1e-5 || norm(K-K_pre) < 1e-5 || iter >= maxIter
    if norm(K-K_pre) < 1e-5 || iter >= maxIter
        terminate = true;
    end
    K_pre = K;

end

% robust regression
maxError = 0.1;
[label, r1, r2] = ssrrr(K,d-1,maxError);
Kt = reduceRankK(K, Vi, D2, Eta, EtaPair, r1, r2, 0.1,maxError);
K = Kt;

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
dimX = nnz(c<0.99)+1;
x = R(1:dimX,:);

x1 = R(1,:);
X1 = hankel(x1(1:d),x1(d:n));
% [Vr, powers] = veronese(X1,2);
[Vr, powers] = veroneseVector(x,d,2);
% Vr2 = reshape(Sh * K(:), p, n-d+1);
Vr2 = K(Vi);
[U,S,V] = svd(Vr);
s = diag(S)
% assert(s(end)<1e-6);
r = U(:,end);

% % my SSC
% Y = hankel(K(1,1:3), K(1,3:end));
% C = ssc(Y);
% W = abs(C) + abs(C.');
% groupSSC1 = SpectralClustering(W, 2);
% % SSC
% OptM = 'L1Perfect'; %OptM can be {'L1Perfect','L1Noise','Lasso','L1ED'}
% CMat = SparseCoefRecovery(Y,0,OptM,0.001);
% CKSym = BuildAdjacency(CMat,0);
% groupSSC2 = SpectralClustering(CKSym,2);


% % vector derivative
% Dpn = zeros(dimX*d,n-d+1);
% for i = 1:n-d+1
%     dv = derivativeVector(r,powers,x(:,i:i+d-1));
%     Dpn(:,i) = dv(:);
% end
% Dpn = normc(Dpn);

%compute only one normal at each point to the subspace containing the point,
%assuming the codim is one.
% Dpn=derivative(c,powers,x);
[Dpn,normDpn] = cnormalize(derivative(r,powers,X1));

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