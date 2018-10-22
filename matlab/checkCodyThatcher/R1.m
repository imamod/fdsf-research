function R_out = R1(p, q, x, n)
    s = 0:n; 
%    temp_nom = 0;
%    temp_denom = 0;
%    for s = 1:n+1
%        temp_nom = temp_nom + p(s)*exp((s-1)*x);
%        temp_denom = temp_denom + q(s)*exp((s-1)*x);
%    end
%    R_out = temp_nom / temp_denom;
    R_out = sum(p.*exp(s*x))/sum(q.*exp(s*x));
end