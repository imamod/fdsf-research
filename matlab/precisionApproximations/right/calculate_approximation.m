function F = calculate_approximation(y, N)

for j = 1:length(y)
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*y(j).^(2*(n-1));
    end
    for m = 1:N
        S2 = S2 + b(m).*y(j).^(2*(m-1));
    end
    F_base(j) = (y(j).^(2*N+2) + S1)/(y(j).^(2*N) + S2);
    delta_base(j) = F_base(j)/z(j)-1;
end