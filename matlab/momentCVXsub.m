function [m, rankCondition, sparseCondition, W1, Ws] = momentCVXsub(opt, var, W1, Ws)

nSys = opt.nSys;
D2 = opt.D2;
D1 = opt.D1;
dataDim = opt.dataDim;
n = opt.n;
d = opt.d;
Eta = opt.Eta;
EtaPair = opt.EtaPair;
lambda1 = opt.lambda1;
epsilon = opt.epsilon;
maxIter = opt.maxIter;
c1 = opt.c1;
c2 = opt.c2;
L = opt.L;
Mi = opt.Mi;
Qi = opt.Qi;
Ki = opt.Ki;
KRi = opt.KRi;
RRi = opt.RRi;
SSi = opt.SSi;
% K2i = opt.K2i;
R2i = opt.R2i;
R1i = opt.R1i;
% Pi = opt.Pi;
Si = opt.Si;

h = opt.h;

eps = 1e-10;

if nargin < 3
    W1 = eye(h);    % weights for moment matrix
    Ws = ones(size(Si));    % weights for sparse indicator variables
end
m = zeros(L, 1);
m_pre = ones(L, 1);
sigma = inf;
terminate = false;
iter = 1;    
while ~terminate
    cvx_clear;
    cvx_begin quiet
    cvx_solver mosek
        variables m(L, 1);
        K = m(Ki);
        for i = 1:n
            for j = 1:n
                if Eta(i,j)==1 || EtaPair(i,j)==1
                    if epsilon==0
                        K(i,i)+K(j,j)-2*K(i,j) >= c1^2*D2(i,j);
                        K(i,i)+K(j,j)-2*K(i,j) <= c2^2*D2(i,j);
                    else
                        K(i,i)+K(j,j)-2*K(i,j) >= c1^2*( D2(i,j) - 4*epsilon*D1(i,j) );
                        K(i,i)+K(j,j)-2*K(i,j) <= c2^2*(D2(i,j) + 4*dataDim*epsilon^2 + 4*epsilon*D1(i,j));
                    end
                end
            end
        end
        K == semidefinite(n);
%         sum(K(:)) == 0;
        m(Mi(1,1)) == 1;
        m(Mi) == semidefinite(size(Mi,1));
%         for i = 1:size(Pi, 3)
%             m(Pi(:,:,i)) == semidefinite(size(Pi,1));
%         end
        
        sum(m(KRi), 2) <= eps;
        sum(m(KRi), 2) >= -eps;
        
        m(RRi) * [-1; 1; 1] <= eps;
        m(RRi) * [-1; 1; 1] >= -eps;
        m(1+var.rBase+1) >= m(1+var.rBase+1+var.rDim);
        sum(m(R2i), 2) == 1;
        m(R1i) >= eps;
        
        m(SSi) * [-1; 1] == 0;
        
        Si2 = reshape(SSi(:,1), nSys, n-d+1);
        sum(m(Si2), 1) == 1;
%         m(SAi) * [-1; 1; 1] <= eps;
%         m(SAi) * [-1; 1; 1] >= -eps;

%         sum(m(R2i)) <= 3;
        
        objMomentRank = trace(W1*m(Qi));
        objSparseS = sum(sum(Ws.*m(Si)));
        obj = objMomentRank + lambda1 * objSparseS;
%         obj = objMomentRank;
%         obj = sum(m(K2i)) + lambda1 * (trace(W1*P) + trace(W2*Q));
        minimize(obj);
    cvx_end
    
    % reweight moment matrix
    sm = svd(m(Qi));
    sigma = min([sigma, sm(2)]);
    W1 = inv(m(Qi) + sigma * eye(h));
    W1 = W1 / norm(W1,'fro');
    % reweight indicator variables to make it sparse
    Ws = 1 ./ (m(Si) + 1e-5);
    Ws = Ws./repmat(sum(Ws),2,1);
    
    iter
    lambda1
    
    S = m(Si);
    S
%     s = svd(m(Mi));
    s = svd(m(Qi));
    s(1:5)'
    
    norm(m-m_pre)/L
    
    rankCondition = (s(2) < 1e-4*s(1));
    sparseCondition = S(1,:)*S(2,:).'<1e-3;
    if lambda1 == 0
        terminate = (rankCondition || norm(m-m_pre)/L<1e-6 || iter >= maxIter);
    else
        terminate = ( (rankCondition && sparseCondition) || norm(m-m_pre)/L<1e-5 || iter >= maxIter );
    end
    m_pre = m;
    iter = iter + 1;
end

end