function Kt = reduceRankK(K, Vi, D2, Eta, EtaPair, r1, r2, errorPercent,maxError)

assert(length(r1)==3 && length(r2)==3);

epsilon = maxError;
maxIter = 100;
errorBound = errorPercent * norm(K,'fro');

n = size(K, 1);

r = zeros(6, 1);
r(1) = r1(1) * r2(1);
r(2) = r1(1) * r2(2) + r1(2) * r2(1);
r(3) = r1(1) * r2(3) + r1(3) * r2(1);
r(4) = r1(2) * r2(2);
r(5) = r1(2) * r2(3) + r1(3) * r2(2);
r(6) = r1(3) * r2(3);


W = eye(n);
Kt_pre = zeros(n, n);
iter = 0;
terminate = false;
while ~terminate
    cvx_begin
    cvx_solver mosek
        variables Kt(n, n);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
                    Kt(i,i)+Kt(j,j)-2*Kt(i,j) <= 25*D2(i,j);
                    Kt(i,i)+Kt(j,j)-2*Kt(i,j) >= 16*D2(i,j);
                end
            end
        end
        Kt == semidefinite(n);
        Kt(Vi) * r <= epsilon;
        Kt(Vi) * r >= -epsilon;
        norm(K-Kt,'fro') <= errorBound;
        obj = trace(W * Kt);
        minimize(obj);
    cvx_end
    
    skt = svd(Kt);
    W = inv(Kt + skt(2)*eye(n));
    W = W / norm(W, 'fro');
    iter = iter + 1;
    
    if norm(Kt-Kt_pre,'fro')<1e-3 || iter==maxIter
        terminate = true;
    end
    Kt_pre = Kt;
end

if iter==maxIter
    fprintf('Problem did not converge.\n');
end

end