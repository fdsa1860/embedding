function Pi = getPSDInd11(Mi, d, n, var)
% get indices to smaller positive definite matrices

nk = var.nk;
nrd = var.nrd;
rDim = var.rDim;
nSys = var.nSys;
kInd = var.kInd;
rInd = var.rInd;
rdInd = var.rdInd;
sInd = var.sInd;

t = n - d + 1;
sharedInd = [1, kInd, rdInd];
m = 1 + nk + nrd + rDim + nSys;
Pi = zeros(m, m, t);
for i = 1:n-d+1
    currInd = [sharedInd, rInd((i-1)*rDim+1:i*rDim), sInd((i-1)*nSys+1:i*nSys)];
    Pi(:,:,i) = Mi(currInd, currInd);
end

end