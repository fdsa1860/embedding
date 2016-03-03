function [X, group] = veroneseSDEcvx11(data, nn, lambda1, lambda2)
% solve 2 systems switch system with moment
% use kernel and regressor and indicator variables
% r = s1 * r1 + s2 * r2;

D = pdist2(data',data');
D = D.^2;
n = size(data, 2);
Eta = getNNmap(D, nn);
EtaPair = (Eta'*Eta > 0);

epsilon = 1e-10;
sysOrd = 2;
d = sysOrd + 1;
nSys = 2;

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
L = max(Indx);
[basis,~] = momentPowers(0,var.nTot,mord);
Mi = getMomInd(Dict,basis,0,Indx,0);
% get indices for regressors
qInd = [1, 1+var.kInd, 1+var.rdInd];
Qi = Mi(qInd, qInd);
% Qi = Mi;
h = size(Qi,1);

% get Kernel indices
Ki = getKernelInd(n);

% Vi = getStructuredVeroneseMap6(Dict, d, n, 2);

KRi = getKRInd11(Dict, d, n, var);
RRi = getRRInd11(Dict, d, n, nSys, var);
SSi = getSSInd11(Dict, var);
% K2i = getK2Ind11(Dict, var);
R2i = getR2Ind11(Dict, var);
R1i = getR1Ind11(Dict, var);

Si = reshape(SSi(:,2), nSys, n-d+1);
% SAi = getSAInd11(Mi, Si);

Ws = ones(size(Si));    % weights for sparse indicator variables

W1 = eye(h);    % weights for moment matrix
m = zeros(L, 1);
m_pre = ones(L, 1);
terminate = false;
maxIter = 1000;
iter = 0;
while ~terminate
    cvx_begin quiet
    cvx_solver mosek
        variables m(L, 1);
        K = m(Ki);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
%                     K(i,i)+K(j,j)-2*K(i,j) >= D(i,j)-epsilon;
%                     K(i,i)+K(j,j)-2*K(i,j) <= D(i,j)+epsilon;
                    K(i,i)+K(j,j)-2*K(i,j) <= 25*D(i,j);
                    K(i,i)+K(j,j)-2*K(i,j) >= 16*D(i,j);
                end
            end
        end
        K == semidefinite(n);
%         sum(K(:)) == 0;
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));
        
        sum(m(KRi), 2) <= epsilon;
        sum(m(KRi), 2) >= -epsilon;
        
        m(RRi) * [-1; 1; 1] <= epsilon;
        m(RRi) * [-1; 1; 1] >= -epsilon;
        m(1+var.rBase+1) >= m(1+var.rBase+1+var.rDim);
        sum(m(R2i), 2) == 1;
        m(R1i) >= epsilon;
        
        m(SSi) * [-1; 1] == 0;
        
        Si2 = reshape(SSi(:,1), nSys, n-d+1);
        sum(m(Si2), 1) == 1;
%         m(SAi) * [-1; 1; 1] <= epsilon;
%         m(SAi) * [-1; 1; 1] >= -epsilon;

%         sum(m(R2i)) <= 3;
        
        objSparseS = sum(sum(Ws.*m(Si)));
        objMomentRank = trace(W1*m(Qi));
        obj = objSparseS + lambda1 * objMomentRank;
%         obj = objMomentRank;
%         obj = sum(m(K2i)) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    % reweight moment matrix
    sm = svd(m(Qi));
    W1 = inv(m(Qi) + sm(2)*eye(h));
    W1 = W1 / norm(W1,'fro');
    % reweight indicator variables to make it sparse
    Ws = 1 ./ (m(Si) + epsilon);
    Ws = Ws./repmat(sum(Ws),2,1);
    
    S = m(Si);
    S
%     s = svd(m(Mi));
    s = svd(m(Qi));
    s(1:5)'
    
    norm(m-m_pre)
%     if s(2)<1e-5*s(1) || norm(m-m_pre)<1e-6 || iter > maxIter
    if (s(2)<1e-5*s(1) && S(1,:)*S(2,:).'<1e-3) || norm(m-m_pre)<1e-3 || iter >= maxIter
        terminate = true;
    end
    m_pre = m;
    iter = iter + 1;
end

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
X = R(1:ind,:);

x1 = R(1,:);
X1 = hankel(x1(1:d),x1(d:end));
[Vr, powers] = veronese(X1,2);
[U,S,V] = svd(Vr);
s = diag(S)
% assert(s(end)<1e-6);
r = U(:,end);

%compute only one normal at each point to the subspace containing the point,
%assuming the codim is one.
% Dpn=derivative(c,powers,x);
[Dpn,normDpn] = cnormalize(derivative(r,powers,X1));

method = 'Cos^2';
%compute the similarity matrix
switch method
    case 'Cos'
        affMat=(abs(Dpn'*Dpn));
    case 'Cos^2'
        affMat=(abs(Dpn'*Dpn)).^2;
    case 'Exp_-sin^2'
        affMat=exp((abs(Dpn'*Dpn)).^2-1); % i.e., exp(-sin^2)
end

%segment using spectral clustering
group=SpectralClustering(affMat,2);

end