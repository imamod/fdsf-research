function max_d = approximate(c2, k, x, y, I_precesion)

I_approximate = zeros(length(c2),length(y));
max_d = zeros(1,length(c2));

c1 = 3*(1-2^(-k))/k;
c3 = (k+1)*pi^2/(gamma(k+2)^(6/k));
c4 = 1/(gamma(k+2)^(6/k));
for j = 1:length(c2)
    for i = 1:length(y)
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+c1*y(i)+c2(j)*y(i)^2+c3*y(i)^4 + c4*y(i)^6)^(k/6);
    end
    delta = log(abs(I_approximate(j,:)./I_precesion));

    plot(x, delta, 'k', 'linewidth',2)
    xlabel('x'); ylabel('ln(d)')
    title(['Подгоночный коэффициент: ' sprintf('%f', c2(j))])
    pause(10^-6);
    if (c2(j) > 1.29)
        break;
    end
    max_d(j) = max(abs(delta));
end
end