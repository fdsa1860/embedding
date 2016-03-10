function [x, label, rHat] = veroneseSDEcvx13(data, opt)
% do not use moment, detect switching between two systems by minimizing
% rank of veronese matrix
% including noise in embedded data X

sysOrd = opt.sysOrd;

D = pdist2(data',data');
opt.D2 = D.^2;
opt.n = size(data, 2);
opt.d = sysOrd + 1;
opt.p = opt.d*(opt.d+1)/2;
opt.Vi = getVeroneseMap5(opt.n, opt.d);
opt.Eta = getNNmap(D, opt.numNeighbors);
opt.EtaPair = (opt.Eta' * opt.Eta > 0);

currKn = [];
opt.lambda1 = opt.lambda1Init;
currLambda1 = opt.lambda1;
terminate2 = false;
while ~terminate2
    [Kn, VeroneseCondition] = veroneseSDEsubNoisy(opt);
    if VeroneseCondition && opt.lambda1 <= 1e6
        currKn = Kn;
        currLambda1 = opt.lambda1;
        opt.lambda1 = opt.lambda1 * opt.lambda1Rate;
    else
        terminate2 = true;
        if ~isempty(currKn)
            Kn = currKn;
            opt.lambda1 = currLambda1;
        end
    end
end

K = Kn(1:opt.n, 1:opt.n);
[x, group, rHat] = gpcaClustering(K, opt);

% assign the first label to initial data
label = [group(1) * ones(1, sysOrd), group.'];
if label(1)~=1
    label = mod(label, opt.nSys) + 1;
end

end