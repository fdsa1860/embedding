% function Y = mySDEcvx(X, k)
function Y = mySDEcvx(D, k)

% d = size(X, 1);
% n = size(X, 2);
% G = X' * X;
n = size(D, 1);
Eta = getNNmap(D, k);
EtaPair = Eta' * Eta;

cvx_begin
cvx_solver sedumi
    variables K(n, n);
    K == semidefinite(n);
    obj = - trace(K);
    sum(K(:)) == 0;
    for i = 1:n
        for j = 1:n
            if Eta(i,j)==1 || EtaPair(i,j)==1
%                 K(i,i)+K(j,j)-2*K(i,j) == G(i,i)+G(j,j)-2*G(i,j);
                K(i,i)+K(j,j)-2*K(i,j) == D(i,j);
            end
        end
    end
    minimize(obj);
cvx_end
[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
Y = R(1:ind,:);

end