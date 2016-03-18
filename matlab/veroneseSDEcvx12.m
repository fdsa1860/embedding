function [x, label, rHat] = veroneseSDEcvx12(data, opt)
% do not use moment, detect switching between two systems by minimizing
% rank of veronese matrix
% automatically select lambda by line search

sysOrd = opt.sysOrd;
nSys = opt.nSys;

D = pdist2(data',data');
opt.D2 = D.^2;
opt.D1 = pdist2(data',data','cityblock');
opt.dataDim = size(data, 1);
opt.n = size(data, 2);
opt.d = sysOrd + 1;
opt.p = opt.d*(opt.d+1)/2;
opt.Vi = getVeroneseMap5(opt.n, opt.d);
opt.Eta = getNNmap(D, opt.numNeighbors);
opt.EtaPair = (opt.Eta' * opt.Eta > 0);

% warm start
opt.lambda1 = 0;
[K, VeroneseCondition, W1,W2] = veroneseSDEsub(opt);
if ~VeroneseCondition
    error('Veronese map not rank deficient.\n');
end
WK = eye(opt.n);

currK = [];
opt.lambda1 = opt.lambda1Init;
currLambda1 = opt.lambda1;
terminate2 = false;
while ~terminate2
%     [K, VeroneseCondition,W1,W2,WK] = veroneseSDEsub(opt,W1,W2,WK);
    [K, VeroneseCondition] = veroneseSDEsub(opt,W1,W2,WK);
%     [K, VeroneseCondition] = veroneseSDEsub(opt);
    if VeroneseCondition && opt.lambda1 <= 1e6
        s = svd(K);
        s(2)/s(1)
        if s(2)/s(1) < 1e-5 % K is rank 1
            terminate2 = true;
        else
            currK = K;
            currLambda1 = opt.lambda1;
            opt.lambda1 = opt.lambda1 * opt.lambda1Rate;
        end
    else
        terminate2 = true;
        if ~isempty(currK)
            K = currK;
            opt.lambda1 = currLambda1;
        end
    end
end

[x, group, rHat] = gpcaClustering(K, opt);

% % robust regression
% maxError = 0.005;
% [group, r1, r2] = ssrrr(K, sysOrd, maxError);
% rHat = [r1 r2];
% x = [];

% % my SSC
% Y = hankel(K(1,1:3), K(1,3:end));
% C = ssc(Y);
% W = abs(C) + abs(C.');
% groupSSC1 = SpectralClustering(W, 2);
% % SSC
% OptM = 'L1Perfect'; %OptM can be {'L1Perfect','L1Noise','Lasso','L1ED'}
% CMat = SparseCoefRecovery(Y,0,OptM,0.001);
% CKSym = BuildAdjacency(CMat,0);
% groupSSC2 = SpectralClustering(CKSym,2);

% assign the first label to initial data
label = [group(1) * ones(1, sysOrd), group.'];
zI = (label==0);
if label(1)~=1
    label = mod(label, nSys) + 1;
end
label(zI) = 0;

end