function Vi = getVeroneseMap8(n, d)
% for nSys = 2 and model noise, get veronese map indices

% p = 2*d*(2*d+1)/2;
p = d*(d+1)/2;
Vi = zeros(n-d+1, p);
for i = 1:n-d+1
    A = tril(true(2*n, 2*n));
    B = false(2*n, 2*n);
    B(i:i+d-1, i:i+d-1) = true;
    C = A & B;
    ind1 = find(C);
%     B = false(2*n, 2*n);
%     B(i+n:i+n+d-1, i+n:i+n+d-1) = true;
%     C = A & B;
%     ind2 = find(C);
%     B = false(2*n, 2*n);
%     B(i+n:i+n+d-1, i:i+d-1) = true;
%     C = A & B;
%     ind3 = find(C);
%     ind = [ind1; ind2; ind3];
    ind = ind1;
    Vi(i,:) = ind.';
end

end