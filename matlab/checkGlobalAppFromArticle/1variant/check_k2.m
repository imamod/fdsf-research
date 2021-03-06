%% �������� ���������� ������������� �� ������
clc
clear all
close all

format long

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
a1 = 0.7:0.01:1;
I_precesion_minus = zeros(length(a1),length(y));
I_precesion = zeros(length(a1),length(y));
I_approximate = zeros(length(a1),length(y));
I_approximate_minus = zeros(length(a1),length(y));
max_d = zeros(1,length(a1));
for j = 1:length(a1)
    for i = 1:length(y)
        I_precesion_minus(j,i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
        I_precesion(j,i) = x(i)^3/3 + x(i)*pi^2/3 + I_precesion_minus(j,i);
        I_approximate(j,i) = gamma(k+1)*y(i).*(1+a1(j)*y(i)+y(i)^3/(gamma(k+2))^(3/k))^(k/3);
        I_approximate_minus(j,i) = gamma(k+1)*y_minus(i).*(1+a1(j)*y_minus(i)+y_minus(i)^3/(gamma(k+2))^(3/k))^(k/3);
    end
    delta = log10(abs(I_approximate(j,:)./I_precesion(j,:)));
    delta_minus = log10(abs(I_approximate_minus(j,:)./I_precesion_minus(j,:)));
    
    d = [delta_minus, delta];
    plot([x_minus, x], d, 'k', 'linewidth', 2)
    axis([-10 30 -0.015 0.015]),grid on
    
%     title(['����������� �����������: ' sprintf('%f', a1(j))])
    pause(10^-6);
    if (a1(j) > 0.76) % a1 = 1.05, mb 1.06
        break;
    end
   
    max_d(j) = max(abs(d));
end
hold on
xlabel('x'); ylabel('lg(I_�_�_�_�_�/I_�_�_�_�)')
plot([-10,0;30,0],[0,-0.015;0,0.015],'k','linewidth',1.4);
disp(max_d');