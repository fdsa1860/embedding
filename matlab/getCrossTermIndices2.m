function ind = getCrossTermIndices2(Dict,d,i,j)

s = sum(Dict, 2);
ind = zeros(d, 1);
for k = 1:d
    if i~=j
        ind(k) = find(Dict(:,i+k-1)==1 & Dict(:,j+k-1)==1 & s==2);
    else
        ind(k) = find(Dict(:,i+k-1)==2 & s==2);
    end
end

end