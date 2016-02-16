function Ki = getKernelInd(n)

Ki = zeros(n,n);
count = 1;
for i = 1:n
    for j = 1:n
        if j >= i
            Ki(i,j) = count;
            count = count + 1;
        else
            Ki(i,j) = Ki(j,i);
        end
    end
end
Ki = Ki + 1;

end