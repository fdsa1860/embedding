function [x,group] = veroneseSDEcvx5(D, k, lambda1, lambda2)
% do not use moment, detect switching between two systems by minimizing
% rank of veronese matrix


D = D.^2;
% d = size(X, 1);
% n = size(X, 2);
% G = X' * X;
n = size(D, 1);
Eta = getNNmap(D, k);
EtaPair = (Eta'*Eta > 0);

epsilon = 1e-6;
d = 3;
p = d*(d+1)/2;

% Sh = getStructuredVeroneseMap5(n, d, 2);
Vi = getVeroneseMap5(n, d);

W1 = eye(n-d+1);
W2 = eye(p);
K = zeros(n, n);
K_pre = ones(n,n);
terminate = false;
while ~terminate
    K_pre = K;
    cvx_begin sdp
    cvx_solver sedumi
        variables K(n, n);
        variables P(n-d+1, n-d+1) Q(p, p);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
%                     K(i,i)+K(j,j)-2*K(i,j) == D(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) >= D(i,j)-epsilon;
                    K(i,i)+K(j,j)-2*K(i,j) <= D(i,j)+epsilon;
                end
            end
        end
        K == semidefinite(n);
        sum(K(:)) == 0;
%         Vr = reshape(Sh * K(:), p, n-d+1).';
        Vr = K(Vi);

        [P Vr; Vr' Q] == semidefinite(p+n-d+1);
%         obj =  - trace(K) + lambda1 * norm_nuc(m(Mi)) +  lambda2 * norm_nuc(m(Vi));
%             obj = - lambda2 * trace(K) ;
        obj = - trace(K) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    s = svd(Vr);
    if s(end) < 1e-5 * s(end-1) || norm(K-K_pre) < 1e-5
        terminate = true;
    end
    
    sp = svd(P);
    sq = svd(Q);
    W1 = inv(P + sp(6)*eye(n-d+1));
    W2 = inv(Q + sq(6)*eye(p));
    scale = norm(blkdiag(W1,W2),2);
    W1 = W1/scale;
    W2 = W2/scale;
end

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
x = R(1:ind,:);

X = hankel(x(1:d),x(d:n));
[Vr, powers] = veronese(X,2);
% Vr2 = reshape(Sh * K(:), p, n-d+1);
Vr2 = K(Vi);
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