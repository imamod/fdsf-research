%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format long

k = 4; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.0560148791230902149024568 ...
     0.0351117957891800867706741 ...
     0.0021834386943672331415760 ...
     0.0002464861525522946634693 ...
     0.0000092228177886669241259];
b = [1 ...
     -0.0611726208769112866900252 ...
      0.0279968542816146833953639 ...
     -0.0007512148294307540141223 ...
      0.0000860680747142919882956];

x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));
a1 = 0.25:0.01:5;
I_precesion_minus = zeros(length(a1),length(y));
I_precesion = zeros(length(a1),length(y));
I_approximate = zeros(length(a1),length(y));
I_approximate_minus = zeros(length(a1),length(y));

max_d = zeros(1,length(a1));
for j = 1:length(a1)
    for i = 1:length(y)
        I_precesion_minus(j,i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
        I_precesion(j,i) = x(i)^5/5 + 2*(x(i)^3)*pi^2/3 + 7*x(i)*pi^4/15 + I_precesion_minus(j,i);
        
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+a1(j)*y(i)+y(i)^3/(gamma(k+2))^(3/k))^(k/3);
        I_approximate_minus(j,i) = gamma(k+1)*y_minus(i).*(1+a1(j)*y_minus(i)+y_minus(i)^3/(gamma(k+2))^(3/k))^(k/3);
    end
    delta = log10(abs(I_approximate(j,:)./I_precesion(j,:)));
    delta_minus = log10(abs(I_approximate_minus(j,:)./I_precesion_minus(j,:)));
    
    d = [delta_minus, delta];
    plot([x_minus, x], [delta_minus, delta], 'k', 'linewidth',2)
    axis([-10 20 -0.03 0.03]),grid on
    xlabel('x'); ylabel('lg(d)')
    title(['Подгоночный коэффициент: ' sprintf('%f', a1(j))])
    pause(10^-6);
    if (a1(j) > 0.48) % a1 = 0.56
        break;
    end
    max_d(j) = max(abs(d));
end
disp(max_d');