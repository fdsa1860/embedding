function [x, group, Dpn] = gpcaClustering(K, opt)

Vi = opt.Vi;
d = opt.d;
n = opt.n;

[U,S,V] = svd(K);
R = S.^0.5 * V';
s = diag(S);
c = cumsum(s)/sum(s);
ind = nnz(c<0.99)+1;
x = R(1:ind,:);

x1 = R(1,:);
X1 = hankel(x1(1:d),x1(d:n));
% [Vr, powers] = veronese(X1,2);
[Vr, powers] = veroneseVector(x,d,2);
Vr2 = K(Vi);
[U,S,V] = svd(Vr);
s = diag(S)
% assert(s(end)<1e-6);
r = U(:,end);

% % vector derivative
% Dpn = zeros(dimX*d,n-d+1);
% for i = 1:n-d+1
%     dv = derivativeVector(r,powers,x(:,i:i+d-1));
%     Dpn(:,i) = dv(:);
% end
% Dpn = normc(Dpn);

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