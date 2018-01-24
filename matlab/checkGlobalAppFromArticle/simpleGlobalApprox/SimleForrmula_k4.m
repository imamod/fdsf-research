%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format shortE

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

I_1_0 = pi^2/12;
I_3_0 = 7*pi^4/120;

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^5 / 5 + 2*x(i)*x(i)*x(i)*pi^2/3 + 7*x(i)*pi^4/15 + I_precesion_minus(i);
    
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