function [V, powers] = veroneseVector(X, d, nSys)
% vector version of veronese map
% Input:
% X: m x n matrix, input data
% d: system order + 1
% nSys: number of systems
% Output:
% V: matrix, veronese map output

[m, n] = size(X);
powers = exponent(2,d);
mV = size(powers, 1);
nV = n - d + 1;
V = zeros(mV, nV);
for i = 1:nV
    for j = 1:mV
        pwr = bsxfun(@times,powers(j,:),ones(m,1));
        t1 = X(:,i:i+d-1).^pwr;
        t2 = prod(t1, 2);
        V(j,i) = sum(t2);
    end
end

end