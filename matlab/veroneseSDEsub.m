function [K, VeroneseCondition,W1,W2,WK] = veroneseSDEsub(opt,W1,W2,WK)

D2 = opt.D2;
D1 = opt.D1;
dataDim = opt.dataDim;
n = opt.n;
d = opt.d;
p = opt.p;
Eta = opt.Eta;
EtaPair = opt.EtaPair;
Vi = opt.Vi;
lambda1 = opt.lambda1;
epsilon = opt.epsilon;
maxIter = opt.maxIter;
c1 = opt.c1;
c2 = opt.c2;

eps = 1e-6;

if nargin < 2
    W1 = eye(n-d+1);
    W2 = eye(p);
    WK = eye(n);
end
K = zeros(n, n);
K_pre = zeros(n,n);
iter = 1;
sigma = inf;
terminate = false;
while ~terminate
    cvx_begin quiet
    cvx_solver mosek
    variables K(n, n);
    variables P(n-d+1, n-d+1) Q(p, p);
    for i = 1:n
        for j = 1:n
            if Eta(i,j)==1 || EtaPair(i,j)==1
                if epsilon==0
%                     K(i,i)+K(j,j)-2*K(i,j) == c1^2*D2(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) >= c1^2*D2(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) <= c2^2*D2(i,j);
                else
                    K(i,i)+K(j,j)-2*K(i,j) >= c1^2*( D2(i,j) - 4*epsilon*D1(i,j) );
                    K(i,i)+K(j,j)-2*K(i,j) <= c2^2*( D2(i,j) + 4*dataDim*epsilon^2 + 4*epsilon*D1(i,j) );
                end
            end
        end
    end
    K == semidefinite(n);
    %         sum(K(:)) == 0;
    Vr = K(Vi);
    
    [P Vr; Vr' Q] == semidefinite(p+n-d+1);
    
    %         obj = (trace(W1*P) + trace(W2*Q)) ;
    obj = (trace(W1*P) + trace(W2*Q)) + lambda1 * trace(WK*K);
    minimize(obj);
    cvx_end
    
    lambda1
    iter
    s = svd(Vr);
    s'
    fprintf('s(end)/sum(s) = %f\n',s(end)/sum(s));
    sk = svd(K);
    sk'
    fprintf('norm(K-K_pre)/n^2 = %f\n',norm(K-K_pre)/n^2);
    
    sigma = min([sigma, s(p)]);
    sigma = max([sigma, eps]);
    W1 = inv(P + sigma * eye(n-d+1));
    W2 = inv(Q + sigma * eye(p));
    scale = norm(blkdiag(W1,W2),2);
    W1 = W1/scale;
    W2 = W2/scale;
    
    WK = inv(K + sk(2)*eye(n));
    WK = WK / norm(WK);
    
    VeroneseCondition = (s(end)/sum(s) < 1e-4);
    if VeroneseCondition || norm(K-K_pre)/n^2 < 1e-5 || iter >= maxIter
        terminate = true;
    end
    K_pre = K;
    iter = iter + 1;
end

end