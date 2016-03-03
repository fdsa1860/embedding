function Ki = getKernelInd4(Dict, n)

Ki = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i==j
        Ki(i, j) = find(sum(Dict, 2)==2 & Dict(:, i)==2);
        else
        Ki(i, j) = find(sum(Dict, 2)==2 & Dict(:, i)==1 & Dict(:, j)==1);
        end
    end
end

end