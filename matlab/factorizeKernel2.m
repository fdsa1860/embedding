function [Rn,R] = factorizeKernel2(K, sysOrd)

n = size(K, 1);
[U,S,V] = svd(K);
N = U(:,sysOrd+1:n);
cvx_begin
cvx_solver mosek
variables P(n-sysOrd, n-sysOrd)
variables Q(n, n-sysOrd)
variables a1 a2 a3 a4
% N * P == Q;
for i = 1:size(Q, 1)
    for j = 1:size(Q, 2)
        if i==j
            Q(i,j) == 1;
        elseif i < j
            Q(i,j) == 0;
        elseif i > j+sysOrd
            Q(i,j) == 0;
%         elseif i == j+1 
%             if j <= 8
%                 Q(i,j) >= a1-1e-3;
%                 Q(i,j) <= a1+1e-3;
%             else
%                 Q(i,j) >= a2-1e-3;
%                 Q(i,j) <= a2+1e-3;
%             end
%         elseif i == j+2
%             if j <= 8
%                 Q(i,j) == a3;
%             else
%                 Q(i,j) == a4;
%             end
        end
    end
end
obj = norm(N * P - Q);
minimize(obj);
cvx_end

R = zeros(sysOrd+1, n-sysOrd);
for i = 1:sysOrd+1
    R(i,:) = diag(Q,-i+1).';
end


Rn = normc(R);

end