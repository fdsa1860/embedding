function Si = getSSInd11(Dict, var)

ns = var.ns; % number of indicator variables
sBase = var.sBase;
Si = zeros(ns, 2);
for i = 1:ns
    Si(i,1) = find(sum(Dict,2)==2 & Dict(:,sBase+i)==2);
    Si(i,2) = find(sum(Dict,2)==1 & Dict(:,sBase+i)==1);
end

end