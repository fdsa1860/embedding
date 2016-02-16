function S = getStructuredVeroneseMap5(n, d, nSys)

p = d*(d+1)/2; % for nSys = 2 only
S = zeros(p*(n-d+1), n^2);
for i = 1:n-d+1
    A = tril(true(n, n));
    B = false(n, n);
    B(i:i+d-1,i:i+d-1) = true;
    C = A & B;
    ind = find(C);
    T = zeros(length(ind), n^2);
    for j = 1:length(ind)
        T(j, ind(j)) = 1;
    end
    S((i-1)*p+1:i*p, :) = T;
end

end