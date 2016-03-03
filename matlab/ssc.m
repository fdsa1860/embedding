function C = ssc(Y)

n = size(Y, 2);
C = zeros(n, n);
for i = 1:n
    y = Y(:, i);
    X = Y(:,[1:i-1 i+1:end]);
    cvx_begin
    cvx_solver mosek
        variables c(n-1)
        y == X * c;
        obj = norm(c, 1);
        minimize(obj);
    cvx_end
    C([1:i-1 i+1:end], i) = c;
end

end