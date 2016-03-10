function [x, label,rHat,rdHat] = veroneseSDEcvx11(data, opt)
% solve 2 systems switch system with moment
% use kernel and regressor and indicator variables
% r = s1 * r1 + s2 * r2;

sysOrd = opt.sysOrd;
nSys = opt.nSys;

d = sysOrd + 1;
n = size(data, 2);
D = pdist2(data',data');

opt.D2 = D.^2;
opt.D1 = pdist2(data',data','cityblock');
opt.dataDim = size(data, 1);
opt.d = d;
opt.n = n;
opt.Eta = getNNmap(D, opt.numNeighbors);
opt.EtaPair = (opt.Eta'*opt.Eta > 0);

var.nSys = nSys;
var.rDim = d;
var.nk = n * (n+1) / 2;  % number of kernel variables
var.nr = var.rDim * (n-d+1);   % number of regressor variables
var.nrd = var.rDim * nSys;     % number of regressor 1 variables
var.ns = nSys * (n-d+1); % number of indicator variables

var.kBase = 0;
var.rBase = var.nk;
var.rdBase = var.rBase + var.nr;
var.sBase = var.rdBase + var.nrd;
var.nTot = var.sBase + var.ns;

var.kInd = (var.kBase+1):(var.kBase+var.nk);
var.rInd = (var.rBase+1):(var.rBase+var.nr);
var.rdInd = (var.rdBase+1):(var.rdBase+var.nrd);
var.sInd = (var.sBase+1):(var.sBase+var.ns);

% build moment matrix
mord = 1;
[Dict,Indx] = momentPowers(0,var.nTot,2*mord);
opt.L = max(Indx);
[basis,~] = momentPowers(0,var.nTot,mord);
Mi = getMomInd(Dict,basis,0,Indx,0);
opt.Mi = Mi;
% get indices for regressors
qInd = [1, 1+var.kInd, 1+var.rdInd];
opt.Qi = Mi(qInd, qInd);
% Qi = Mi;
opt.h = size(opt.Qi, 1);
% get Kernel indices
opt.Ki = getKernelInd(n);

opt.KRi = getKRInd11(Dict, d, n, var);
opt.RRi = getRRInd11(Dict, d, n, nSys, var);
opt.SSi = getSSInd11(Dict, var);
% K2i = getK2Ind11(Dict, var);
opt.R2i = getR2Ind11(Dict, var);
opt.R1i = getR1Ind11(Dict, var);
opt.Pi = getPSDInd11(Mi, d, n, var);

opt.Si = reshape(opt.SSi(:,2), nSys, n-d+1);
% SAi = getSAInd11(Mi, Si);

% % warm start
% opt.lambda1 = 0;
% [m, rankCondition, sparseCondition, W1] = momentCVXsub(opt, var);
% if ~rankCondition
%     error('Moment matrix is not rank 1.\n');
% end
Ws = ones(size(opt.Si));
W1 = eye(opt.h);

% opt.maxIter = 100;
opt.lambda1 = opt.lambda1Init;
terminate2 = false;
while ~terminate2
    Ws = ones(size(opt.Si));
    [m, rankCondition,sparseCondition,W1,Ws] = momentCVXsub(opt,var,W1,Ws);
%     [m, rankCondition, sparseCondition] = momentCVXsub(opt, var, W1, Ws);
%     [K, rankCondition] = momentCVXsub(opt, var);
    if rankCondition && sparseCondition 
        terminate2 = true;
    elseif ~sparseCondition && opt.lambda1 <= 1e6
        opt.lambda1 = opt.lambda1 * opt.lambda1Rate;
    elseif ~rankCondition && opt.lambda1 >= 1e-6
        opt.lambda1 = opt.lambda1 / opt.lambda1Rate;
    end
end

save ../expData/moment_n20_e02_m.mat m;

K = m(opt.Ki);
[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
x = R(1:ind,:);

[~,group] = max(m(opt.Si));

Ri = reshape(var.rInd, var.rDim, []);
Ri = Ri + 1;
rHat = m(Ri);

Rdi = reshape(var.rdInd, var.rDim, var.nSys);
Rdi = Rdi + 1;
rdHat = m(Rdi);

% assign the first label to initial data
label = [group(1) * ones(1, sysOrd), group];
if label(1)~=1
    label = mod(label, opt.nSys) + 1;
end

end