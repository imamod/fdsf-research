%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format long

k = 1; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.2715113138214362780964488 ...
     0.0562661238060587633169245 ...
     0.0067420740469345689743873 ...
     0.0005169505155333205859985 ...
     0.0000194771836765773190602];
b = [1 ...
     0.0215113138214352840651584 ...
     0.0231105175729721417901084 ...
     0.0003669081577365413477999 ...
     0.0000610424408732720110769];

x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^2/2 + pi^2/6 - I_precesion_minus(i);
    
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
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)) + 10^-17);
end

% disp('delta = ');
% disp((delta)');
grid on, hold on
plot(x, delta, '*')
plot(x_minus, delta_minus, '*')