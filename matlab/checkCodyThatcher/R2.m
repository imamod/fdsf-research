function R_out = R2(p, q, x, n)
   s = 0:n;
%    temp = p.*(x.^s)
   R_out = (sum(p.*(x.^s)))/(sum(q.*(x.^s)));
end