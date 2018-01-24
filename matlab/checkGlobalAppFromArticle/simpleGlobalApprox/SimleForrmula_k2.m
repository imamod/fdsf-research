%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format shortE

k = 2; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.2263816364340698560028783 ...
     0.0533684335574798857246766 ...
     0.0062904756340795211604491 ...
     0.0005023228274452983506998 ...
     0.0000189379675088061004880];
 
b = [1 ...
     0.0388816364340691133155655 ...
     0.0243043998742774445085992 ...
     0.0006290985326433190105734 ...
     0.0000657018161945458806177];

x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));

I_1_0 = pi^2/12;

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^3 / 3 + x(i)*pi^2/3 + I_precesion_minus(i);
    
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