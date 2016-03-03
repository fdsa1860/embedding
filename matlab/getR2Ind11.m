function R2i = getR2Ind11(Dict, var)
% find all second order monomials of regressor variables

ind = [var.rInd, var.rdInd];
R2i = zeros(1, length(ind));
for i = 1:length(ind)
    R2i(i) = find(sum(Dict,2)==2 & Dict(:,ind(i))==2);
end
R2i = reshape(R2i, var.rDim, []).';

end