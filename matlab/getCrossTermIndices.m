function ind = getCrossTermIndices(Dict,d,i,j)

s = sum(Dict, 2);
ind = zeros(d, 1);
for k = 1:d
    if i~=j
        ind(k) = find(Dict(:,d*(i-1)+k)==1 & Dict(:,d*(j-1)+k)==1 & s==2);
    else
        ind(k) = find(Dict(:,d*(i-1)+k)==2 & s==2);
    end
end

end