function Y = veroneseSDEglopt(D, k, lambda)

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
    % gloptipoly begin
        mset clear;
        mpol x 10 2;
        K = x*x';

        C = [];
        %         K == semidefinite(n);
        %         C = [C; K >= 0];
        C = [C; sum(K(:)) == 0];
        for i = 1:n
            for j = 1:n
%                 K(i,j) == X(i,2)*X(j,2) + X(i,3)*X(j,3);
                if Eta(i,j)==1 || EtaPair(i,j)==1
                    C = [C; K(i,i)+K(j,j)-2*K(i,j) == D(i,j)];
                end
            end
        end
        g0 = - lambda * trace(K);
        P = msdp(min(g0),C,2);
        [status, obj] = msol(P);

    %     [P G; G' Q] == semidefinite(2*m);
    %     obj = trace(W1*Y) + trace(W2*Z) - lambda * trace(K);
    %     obj = trace(W1*P) + trace(W2*Q);
    % gloptipoly end
    
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