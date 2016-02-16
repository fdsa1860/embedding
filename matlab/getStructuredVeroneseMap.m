function S = getStructuredVeroneseMap(Dict, d, nVar, nSys)

L = size(Dict, 1);
p = d*(d+1)/2; % for nSys = 2 only
S = zeros(p*(nVar-d+1), L);
for i = 1:nVar-d+1
    ind = find(sum(Dict(:,i:i+d-1),2) == nSys & sum(Dict,2) == nSys)';
    T = zeros(length(ind), L);
    for j = 1:length(ind)
        T(j, ind(j)) = 1;
    end
    S((i-1)*p+1:i*p, :) = T;
end

end