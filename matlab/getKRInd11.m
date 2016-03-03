function KRi = getKRInd11(Dict, d, n, var)

Ki = getKernelInd(n) - 1;

KRi = zeros(n*(n-d+1), d);
for i = 1:n-d+1
    for j = 1:n
        for k = 1:d
            KRi((i-1)*n+j, k) = find(sum(Dict,2)==2 & Dict(:,Ki(j,i+k-1))==1 & Dict(:,var.nk+(i-1)*d+k)==1);
%             if k == 1
%                 KRi((i-1)*n+j, k) = find(sum(Dict,2)==1 & Dict(:,Ki(j,i+k-1))==1);
%             else
%                 KRi((i-1)*n+j, k) = find(sum(Dict,2)==2 & Dict(:,Ki(j,i+k-1))==1 & Dict(:,var.nk+(i-1)*d+k-1)==1);
%             end
        end
    end
end

end