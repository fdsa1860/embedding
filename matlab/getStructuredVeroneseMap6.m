function Vi = getStructuredVeroneseMap6(Dict, d, n, nSys)

q = n*(n+1)/2;
p = d*(d+1)/2;
Ki = getKernelInd(n);
Ki = Ki - 1;

Vi = zeros(n-d+1, p);
for k = 1:n-d+1
    cnt = 1;
    for i = k:k+d-1
        for j = i:k+d-1
            Vi(k,cnt) = find(Dict(:,Ki(i,j))==1 & Dict(:,q+cnt)==1 & sum(Dict,2)==nSys);
            cnt = cnt + 1;
        end
    end
end

end