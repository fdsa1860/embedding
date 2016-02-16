function S = getStructuredVeroneseMap4(Dict, d, nVar, nSys)
% get monomials of veronese map times coefficients

L = size(Dict, 1);
p = d*(d+1)/2; % for nSys = 2 only
S = zeros(p*(nVar-d+1), L);
for i = 1:nVar-d+1
    Index = find(sum(Dict(:,i:i+d-1),2) == nSys & sum(Dict(:,1:nVar),2) == nSys & sum(Dict,2) == nSys+1)';
    Index = reshape(Index,[p, p]);
    ind = diag(Index);
    T = zeros(length(ind), L);
    for j = 1:length(ind)
        T(j, ind(j)) = 1;
    end
    S((i-1)*p+1:i*p, :) = T;
end

end