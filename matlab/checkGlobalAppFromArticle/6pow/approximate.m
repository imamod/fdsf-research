function max_d = approximate(a, b, c2, k, P, x, y, x_minus, y_minus, N)

n = 0:N+1; m = 0:N;

I_precesion_minus = zeros(length(c2),length(y));
I_precesion = zeros(length(c2),length(y));
I_approximate = zeros(length(c2),length(y));
I_approximate_minus = zeros(length(c2),length(y));
max_d = zeros(1,length(c2));

c1 = 3*(1-2^(-k))/k;
c3 = (k+1)*pi^2/(gamma(k+2)^(6/k));
c4 = 1/(gamma(k+2)^(6/k));
for j = 1:length(c2)
    for i = 1:length(y)
        I_precesion_minus(j,i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
        if (rem(k,2) == 0)
            I_precesion(j,i) = P(i) + I_precesion_minus(j,i);
        else
            I_precesion(j,i) = P(i) - I_precesion_minus(j,i);
        end
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+c1*y(i)+c2(j)*y(i)^2+c3*y(i)^4 + c4*y(i)^6)^(k/6);
        I_approximate_minus(j,i) = gamma(k+1)*y_minus(i).*(1+c1*y_minus(i)+c2(j)*y_minus(i)^2+c3*y_minus(i)^4 + c4*y_minus(i)^6)^(k/6);
    end
    delta = log10(abs(I_approximate(j,:)./I_precesion(j,:)));
    delta_minus = log10(abs(I_approximate_minus(j,:)./I_precesion_minus(j,:)));
    
    d = [delta_minus, delta];
    plot([x_minus, x], [delta_minus, delta], 'k', 'linewidth',2)
%     axis([-10 30 -0.01 0.01]),grid on
    xlabel('x'); ylabel('lg(d)')
    title(['Подгоночный коэффициент: ' sprintf('%f', c2)])
    pause(10^-6);
%     if (c2(j) > 1)
%         break;
%     end
    max_d(j) = max(abs(d));
end
end