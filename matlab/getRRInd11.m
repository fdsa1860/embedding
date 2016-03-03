function RRi = getRRInd11(Dict, d, n, nSys, var)

m = var.rDim; % regressor dimension
rBase = var.rBase;
rdBase = var.rdBase;
sBase = var.sBase;

RRi = zeros(m*(n-d+1), nSys+1);
for i = 1:n-d+1
    for j = 1:m
        rInd = (i-1)*m+j;
        for k = 1:nSys+1
            if k == 1
                RRi(rInd, k) = find(sum(Dict,2)==1 & Dict(:,rBase+rInd)==1);
            else
                rdInd = (k-2)*m+j;
                sInd = (i-1)*nSys+k-1;
                RRi(rInd, k) = find(sum(Dict,2)==2 & Dict(:,rdBase+rdInd)==1 & Dict(:,sBase+sInd)==1);
            end
        end
    end
end

end