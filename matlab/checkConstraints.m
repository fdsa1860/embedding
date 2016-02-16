% checkConstraints

assert(abs(sum(K(:))) < 1e-3);
for i = 1:n
    for j = 1:n
        if Eta(i,j)==1 || EtaPair(i,j)==1
            assert(abs(K(i,i)+K(j,j)-2*K(i,j) - D(i,j))<1e-3);
        end
    end
end

fprintf('Success!\n')