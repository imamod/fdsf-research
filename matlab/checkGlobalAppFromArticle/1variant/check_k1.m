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
a1 = 0.8:0.01:5;
I_precesion_minus = zeros(length(a1),length(y));
I_precesion = zeros(length(a1),length(y));
I_approximate = zeros(length(a1),length(y));
I_approximate_minus = zeros(length(a1),length(y));
max_d = zeros(1,length(a1));
for j = 1:length(a1)
    for i = 1:length(y)
        I_precesion_minus(j,i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
        I_precesion(j,i) = x(i)^2/2 + pi^2/6 - I_precesion_minus(j,i);
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+a1(j)*y(i)+y(i)^3/(gamma(k+2))^(3/k))^(k/3);
        I_approximate_minus(j,i) = gamma(k+1)*y_minus(i).*(1+a1(j)*y_minus(i)+y_minus(i)^3/(gamma(k+2))^(3/k))^(k/3);
    end
    delta = log10(abs(I_approximate(j,:)./I_precesion(j,:)));
    delta_minus = log10(abs(I_approximate_minus(j,:)./I_precesion_minus(j,:)));
    
    d = [delta_minus, delta];
    plot([x_minus, x], [delta_minus, delta], 'k', 'linewidth',2)
    axis([-10 30 -0.009 0.009]),grid on
    xlabel('x'); ylabel('lg(d)')
    title(['Подгоночный коэффициент: ' sprintf('%f', a1(j))])
    pause(10^-6);
    if (a1(j) > 1)
        break;
    end
    max_d(j) = max(abs(d));
end
hold on
plot([-10,0;30,0],[0,-0.06;0,0.06],'k','linewidth',1.5);
disp(max_d');