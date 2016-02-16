function [Y,group] = veroneseSDEcvx(D, k, lambda1, lambda2)

D = D.^2;
% d = size(X, 1);
% n = size(X, 2);
% G = X' * X;
n = size(D, 1);
Eta = getNNmap(D, k);
EtaPair = Eta' * Eta;

d = 3;
mord = 1;
[Dict,Indx] = momentPowers(0,d*n,2*mord);
L = max(Indx);

% build moment matrix
[basis,~] = momentPowers(0,d*n,mord);

Mi = getMomInd(Dict,basis,0,Indx,0);

Vi = getVeroneseMap(Dict, basis, d, 2);

h = size(Mi,1);

W1 = eye(h);
W2 = eye(h);
m = zeros(L, 1);
m_pre = ones(L, 1);
while norm(m - m_pre) > 1e-3
    m_pre = m;
    cvx_begin sdp
    cvx_solver sedumi
        variables m(L, 1);
        variables K(n, n);
        variables P(h, h) Q(h, h);
        for i = 1:n
            for j = 1:n
                ind = getCrossTermIndices(Dict,d,i,j);
                K(i,j) == sum(m(ind));
                if Eta(i,j)==1 || EtaPair(i,j)==1
                    K(i,i)+K(j,j)-2*K(i,j) == D(i,j);
%                     X(i,1)+X(i,4) + X(j,1)+X(j,4) - sqrt(X(i,2)-X(i,3) - X(j,2)-X(j,3) == D(i,j);
                end
            end
        end
        K == semidefinite(n);
        sum(K(:)) == 0;
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));
        
        H = hankel(1:d,d:n+d-1);
        vec_H = H(:);
        for i = 1:n+d-1
            for j = i+1:n+d-1
                if vec_H(j)==vec_H(i)
                    m(j) == m(i);
                    ind = find(Dict(:,i)>0);
                    for k = 1:length(ind)
                        currRow = Dict(ind(k),:);
                        temp = currRow(j);
                        currRow(j) = currRow(i);
                        currRow(i) = temp;
                        newInd = find(ismember(Dict,currRow,'rows'));
                        m(newInd) == m(ind(k));
                        Dict(newInd,:)
                        Dict(ind(k),:)
                    end
                end
            end
        end

        [P m(Mi); m(Mi)' Q] == semidefinite(2*h);
%         obj =  - trace(K) + lambda1 * norm_nuc(m(Mi)) +  lambda2 * norm_nuc(m(Vi));
%             obj = - lambda2 * trace(K) ;
    %     obj = trace(W1*Y) + trace(W2*Z) - lambda * trace(K);
        obj = - trace(K) + lambda1 * (trace(W1*P) + trace(W2*Q)) + lambda2 * norm_nuc(m(Vi));
        minimize(0);
    cvx_end
    
    sy = svd(P);
    sz = svd(Q);
    W1 = inv(P + sy(2)*eye(h));
    W2 = inv(Q + sz(2)*eye(h));
end

% [U,S,V] = svd(K);
% R = S.^0.5 * V';
% s = diag(S);
% c = cumsum(s)/sum(s);
% ind = nnz(c<0.99)+1;
% Y = R(1:ind,:);

Y = reshape(m(2:d*n+1),d,n);
[Vr, powers] = veronese(Y,2);
[U,S,V] = svd(Vr);
s = diag(S);
assert(s(end)<1e-6);
r = U(:,end);

%compute only one normal at each point to the subspace containing the point,
%assuming the codim is one.
% Dpn=derivative(c,powers,x);
[Dpn,normDpn] = cnormalize(derivative(r,powers,Y));

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