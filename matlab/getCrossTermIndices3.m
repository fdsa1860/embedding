function ind = getCrossTermIndices3(Dict,i,j)

s = sum(Dict, 2);
if i~=j
    ind = find(Dict(:,i)==1 & Dict(:,j)==1 & s==2);
else
    ind = find(Dict(:,i)==2 & s==2);
end

end