function [Kn, VeroneseCondition] = veroneseSDEsubNoisy(opt)

D2 = opt.D2; 
n = opt.n;
d = opt.d;
p = opt.p;
Eta = opt.Eta;
EtaPair = opt.EtaPair;
Vi = opt.Vi;
lambda1 = opt.lambda1;
maxIter = opt.maxIter;
epsilon = opt.epsilon;
c1 = opt.c1;
c2 = opt.c2;

W1 = eye(n-d+1);
W2 = eye(p);
WK = eye(2*n);
Kn = zeros(2*n, 2*n);
Kn_pre = zeros(2*n,2*n);
iter = 1;
terminate = false;
while ~terminate
    cvx_begin quiet
    cvx_solver mosek
    variables Kn(2*n, 2*n);
    variables P(n-d+1, n-d+1) Q(p, p);
    for i = 1:n
        for j = 1:n
            if Eta(i,j)==1 || EtaPair(i,j)==1
                Kn(i,i)+Kn(j,j)+Kn(i+n,i+n)+Kn(j+n,j+n)...
                    +2*Kn(i,i+n)+2*Kn(j,j+n)-2*Kn(i,j)-2*Kn(i+n,j+n)...
                    -2*Kn(i,j+n) - 2*Kn(i+n,j) >= c1^2 * D2(i,j);
                Kn(i,i)+Kn(j,j)+Kn(i+n,i+n)+Kn(j+n,j+n)...
                    +2*Kn(i,i+n)+2*Kn(j,j+n)-2*Kn(i,j)-2*Kn(i+n,j+n)...
                    -2*Kn(i,j+n) - 2*Kn(i+n,j) <= c2^2 * D2(i,j);
            end
        end
    end
    Kn == semidefinite(2*n);
    max(max(Kn(1+n:2*n,1+n:2*n))) <= epsilon^2;
    min(min(Kn(1+n:2*n,1+n:2*n))) >= -epsilon^2;
    %         sum(K(:)) == 0;
    K = Kn(1:n,1:n);
    Vr = K(Vi);
    
    [P Vr; Vr' Q] == semidefinite(p+n-d+1);
    
    %         obj = (trace(W1*P) + trace(W2*Q)) ;
    obj = (trace(W1*P) + trace(W2*Q)) + lambda1 * trace(WK*Kn);
    minimize(obj);
    cvx_end
    
    s = svd(Vr);
    s'
    sk = svd(Kn);
    sk'
    
    sp = svd(P);
    sq = svd(Q);
    W1 = inv(P + sp(5)*eye(n-d+1));
    W2 = inv(Q + sq(5)*eye(p));
    scale = norm(blkdiag(W1,W2),2);
    W1 = W1/scale;
    W2 = W2/scale;
    
    WK = inv(Kn + sk(2)*eye(2*n));
    WK = WK / norm(WK);
    
    VeroneseCondition = (s(end)/sum(s) < 1e-4);
    if VeroneseCondition || norm(Kn-Kn_pre)/n^2 < 1e-5 || iter >= maxIter
        terminate = true;
    end
    Kn_pre = Kn;
    iter = iter + 1;
end

end