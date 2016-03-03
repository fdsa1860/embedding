function K2i = getK2Ind11(Dict, var)
% find all second order monomials of regressor variables

kInd = var.kInd;
K2i = zeros(length(kInd), 1);
for i = 1:length(kInd)
    K2i(i) = find(sum(Dict,2)==2 & Dict(:,kInd(i))==2);
end

end