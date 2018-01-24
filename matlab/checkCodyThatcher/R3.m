function R_out = R3(p, q, x, n)
   s = 0:n;
%    temp = p.*(x.^s)
   R_out = sum(p.*(x.^(-2*s)))/sum(q.*(x.^(-2*s)));
end