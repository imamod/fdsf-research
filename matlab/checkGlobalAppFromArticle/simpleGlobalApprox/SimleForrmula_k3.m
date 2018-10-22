%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format shortE

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

I_1_0 = pi^2/12;
I_3_0 = 7*pi^4/120;

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^4 / 4 + x(i)*x(i)*pi^2/2 + 7*pi^4/60 - I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % Для страховки добавим 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)));
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
    % Для страховки добавим 10^-17
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)));
end

disp('delta = ');
disp((delta)');
grid on, hold on
xlabel('X'); ylabel('lg(d)')
plot(x, delta, '*')
plot(x_minus, delta_minus, '*')