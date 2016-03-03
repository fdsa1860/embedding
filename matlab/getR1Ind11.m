function R1i = getR1Ind11(Dict, var)
% find all second order monomials of regressor variables

ind = [var.rInd, var.rdInd];
Ri = zeros(1, length(ind));
for i = 1:length(ind)
    Ri(i) = find(sum(Dict,2)==1 & Dict(:,ind(i))==1);
end
Ri = reshape(Ri, var.rDim, []).';
R1i = Ri(:, 1);

end