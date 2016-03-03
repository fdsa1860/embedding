function SAi = getSAInd11(Mi, Si)
% get indices of monomials including s times all other variables

nBlock = size(Mi, 1) - 1;
nSys = size(Si, 1);
nData = size(Si, 2);
SAi = zeros(nData*nBlock, 3);
for i = 1:nData
    r = nBlock * (i-1) + 1;
    for j = 1:nSys+1
        if j == 1
            SAi(r:r+nBlock-1, j) = Mi(2:end, 1);
        else
            SAi(r:r+nBlock-1, j) = Mi(2:end, Si(j-1, i));
        end
    end
end

end