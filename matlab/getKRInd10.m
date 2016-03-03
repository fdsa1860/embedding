function KRi = getKRInd10(Dict, d, n)

q = n*(n+1)/2; % number of kernel variables
p = d*(d+1)/2; % number of veronese map variables
n1 = d;
n2 = d;
baseInd = q+p+n1+n2;
Ki = getKernelInd(n) - 1;

KRi = zeros(n*(n-d+1),d);
for i = 1:n-d+1
    for j = 1:n
        for k = 1:d
            KRi((i-1)*n+j, k) = find(sum(Dict,2)==2 & Dict(:,Ki(j,i+k-1))==1 & Dict(:,baseInd+(i-1)*d+k)==1);
        end
    end
end

end