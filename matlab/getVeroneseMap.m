function ind = getVeroneseMap(Dict, basis, d, nSys)

ind = [];
nVar = size(basis, 2);
for i = 1:d:nVar
    ind = [ind; find(sum(Dict(:,i:i+d-1),2) == nSys)'];
end

end