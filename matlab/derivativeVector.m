function [Dp,Dc] = derivativeVector(c,powers,x)
% take derivative wrt vector, generic from Rene Vidal's GPCA code

n = max(powers(:,1));
K = size(powers,2);

x(abs(x) < 1e-10) = 1e-10;

Dc = (c*ones(1,K)).* powers; 
Dc = Dc(:);
Dc = Dc(powers>0);
Dc = reshape(Dc,length(Dc)/K,K);

Dp = x * Dc;


end