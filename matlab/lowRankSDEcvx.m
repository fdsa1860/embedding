function Y = lowRankSDEcvx(D, k, lambda)

D = D.^2;
% d = size(X, 1);
% n = size(X, 2);
% G = X' * X;
n = size(D, 1);
Eta = getNNmap(D, k);
EtaPair = Eta' * Eta;

m = floor(n/2);

% W1 = eye(m);
% W2 = eye(m);
% K = zeros(n, n);
% K_pre = ones(n, n);
% while norm(K - K_pre,'fro') > 1e-3
%     K_pre = K;
    cvx_begin
    cvx_solver sedumi
    variables K(n, n);
%     variables P(m, m) Q(m, m);
    K == semidefinite(n);
    sum(K(:)) == 0;
    for i = 1:n
        for j = 1:n
            if Eta(i,j)==1 || EtaPair(i,j)==1
                K(i,i)+K(j,j)-2*K(i,j) == D(i,j);
            end
        end
    end
    
    G = zeros(m,m);
    for i = 1:n-m
        G = G + K(i:i+m-1,i:i+m-1);
    end
%     [P G; G' Q] == semidefinite(2*m);
        obj = norm_nuc(G) - lambda * trace(K);
%     obj = trace(W1*Y) + trace(W2*Z) - lambda * trace(K);
%     obj = trace(W1*P) + trace(W2*Q);
    minimize(obj);
    cvx_end
    
%     sy = svd(P);
%     sz = svd(Q);
%     W1 = inv(P + sy(2)*eye(m));
%     W2 = inv(Q + sz(2)*eye(m));
% end

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
Y = R(1:ind,:);

end