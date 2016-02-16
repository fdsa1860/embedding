function Vi = getVeroneseMap5(n, d)
% for nSys = 2 and model noise, get veronese map indices

p = d*(d+1)/2; 
Vi = zeros(n-d+1, p);
for i = 1:n-d+1
    A = tril(true(n, n));
    B = false(n, n);
    B(i:i+d-1, i:i+d-1) = true;
    C = A & B;
    ind = find(C);
    Vi(i,:) = ind.';
end

end