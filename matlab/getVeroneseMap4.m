function Vi = getVeroneseMap4(Dict, d, n, nSys)
% for nSys = 2, get veronese map indices

p = d*(d+1)/2; % for nSys = 2 only
Vi = zeros(n-d+1, p);
for i = 1:n-d+1
    Index = find(sum(Dict(:,i:i+d-1),2) == nSys & sum(Dict(:,1:n),2) == nSys & sum(Dict,2) == nSys+1)';
    Index = reshape(Index,[p, p]);
    ind = diag(Index);
    Vi(i,:) = ind.';
end

end