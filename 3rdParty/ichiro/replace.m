function [ ret ] = replace( mat, from, to )
%REPLACE Summary of this function goes here
%   Detailed explanation goes here

m = containers.Map;
for k=1:length(from)
  m(char(from(k)))=to(k);
end
m('1')=1;

  function elem_ret = replace_element(elem)
    [c,t] = coeffs(expand(elem));
    elem_ret = 0;
    for l=1:length(c)
      elem_ret = elem_ret + double(c(l)) * m(char(t(l)));
    end
  end

ret=to(1)*zeros(size(mat));
for k1=1:size(mat,1)
  for k2=1:size(mat,2)
    ret(k1,k2) = replace_element(mat(k1,k2));
  end
end
end

