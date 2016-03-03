function [r1, r2] = factorizeKernel(r, K, sysOrd)

% build moment matrix
mord = 1;
n = size(K, 1);
[Dict,Indx] = momentPowers(0,n*sysOrd,2*mord);
L = max(Indx);
[basis,~] = momentPowers(0,n*sysOrd,mord);
Mi = getMomInd(Dict,basis,0,Indx,0);
Ki = getKernelIndex(Dict, sysOrd);

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
        variables m(L,1);
        variables P(h, h) Q(h, h);
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));


        [P m(Mi); m(Mi)' Q] == semidefinite(2*h);

        obj = (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    s = svd(m(Mi));
    if s(2)<1e-5*s(1) || norm(m-m_pre,inf)<1e-5 || iter > maxIter
        terminate = true;
    end
    iter = iter + 1;
    
    sy = svd(P);
    sz = svd(Q);
    W1 = inv(P + sy(2)*eye(h));
    W2 = inv(Q + sz(2)*eye(h));
    %     scale = norm(blkdiag(W1,W2),2);
    %     W1 = W1/scale;
    %     W2 = W2/scale;
end

svd(m(Mi))
r1 = m(1:3);
r2 = m(4:6);

end