function [ ret ] = monomials( vars, order)
%MONOMIALS Summary of this function goes here
%   Detailed explanation goes here

  function uret=unsorted_monomials(vars, order)
    uret = [];
    if isempty(vars)
      uret =1;
    else
      for k=0:order
        uret=[uret,vars(1)^k*unsorted_monomials(vars(2:end),order-k)];
      end
    end
  end

  ret=unsorted_monomials(vars, order);
  order = log2(subs(ret,vars,2*ones(size(vars))));
  [~,idx] = sort(order);
  ret = ret(idx);
end