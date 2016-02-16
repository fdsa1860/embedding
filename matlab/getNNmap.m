function Eta = getNNmap(G,k)
% Given pairwise distance matrix G, get nearest neighbor map Eta

assert(ismatrix(G));
assert(size(G,2)==size(G,2));

n = size(G,1);
Eta = zeros(n,n);
% [~, ind] = sort(G,'descend');
[~, ind] = sort(G,'ascend');
nnIndex = ind(1:k,:);
for i = 1:n
    Eta(nnIndex(:,i),i) = 1;
end

end