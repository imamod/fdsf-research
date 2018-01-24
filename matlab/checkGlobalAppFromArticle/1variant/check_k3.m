%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format long

k = 3; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.1583482145380455955096383 ...
     0.0460645149909308107878344 ...
     0.0048861379108841469134267 ...
     0.0004336733305971515517559 ...
     0.0000173435613795895152436];
b = [1 ...
     0.0125148812047107612191739 ...
     0.0266693407000929631393759 ...
     0.0003285431094547362504004 ...
     0.0000820910787890062715299];

x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));
a1 = 0.5:0.01:5;
I_precesion_minus = zeros(length(a1),length(y));
I_precesion = zeros(length(a1),length(y));
I_approximate = zeros(length(a1),length(y));
I_approximate_minus = zeros(length(a1),length(y));
max_d = zeros(1,length(a1));
for j = 1:length(a1)
    for i = 1:length(y)
        I_precesion_minus(j,i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
        I_precesion(j,i) = x(i)^4/4 + 6*x(i)*x(i)*pi^2/12 + 7*pi^4/60 - I_precesion_minus(j,i);
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+a1(j)*y(i)+y(i)^3/(gamma(k+2))^(3/k))^(k/3);
        I_approximate_minus(j,i) = gamma(k+1)*y_minus(i).*(1+a1(j)*y_minus(i)+y_minus(i)^3/(gamma(k+2))^(3/k))^(k/3);
    end
    delta = log10(abs(I_approximate(j,:)./I_precesion(j,:)));
    delta_minus = log10(abs(I_approximate_minus(j,:)./I_precesion_minus(j,:)));
    
    d = [delta_minus, delta];
    plot([x_minus, x], [delta_minus, delta], 'k','linewidth' ,2)
    axis([-10 30 -0.025 0.025]),grid on
%     axis([-10 50 -0.3 0.3]),grid on
    xlabel('x'); ylabel('lg(d)')
    title(['Подгоночный коэффициент: ' sprintf('%f', a1(j))])
    pause(10^-6);
    if (a1(j) > 0.59) % a1 = 0.37
        break;
    end
    max_d(j) = max(abs(d));
end
disp(max_d');