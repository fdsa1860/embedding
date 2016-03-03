function [m] = factorizeRegressor8(r, K, sysOrd)

epsilon = 1e-1;
n = size(K, 1);

Vi = getVeroneseMap5(n, sysOrd+1);

cvx_begin sdp
cvx_solver mosek
    variables m(n, n, sysOrd);
    sum(m,3) == K;
    for i = 1:sysOrd
        M = m(:,:,i);
        M(Vi) * r <= epsilon;
        M(Vi) * r >= -epsilon;
%         sum(abs(M(:))) >= 0.0001;
    end
%     m(:) >= 0.01;
    minimize(0);
cvx_end

end