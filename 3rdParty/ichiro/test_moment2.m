syms x1 x2 real

f0 = x1^3 - 2*x2^2
f1 = 1 - x1^2 - x2^2

u1 = monomials([x1,x2],1)'
u2 = monomials([x1,x2],2)'

mono = monomials([x1,x2],4)
mono(1)=[];

cvx_begin sdp
cvx_solver sedumi;
variable m(size(mono));
minimize replace(f0, mono, m)
subject to
replace(u1*u1'*f1, mono, m) >=0 
replace(   u2*u2', mono, m) >=0
cvx_end


f0val = replace(f0,mono,m)
f1val = replace(f1,mono,m)