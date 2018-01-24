%% �������� ���������� ������������� �� ������
clc
clear all
close all

format shortE

k = 1; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.3126028287472988...
     0.0673008212829461...
     0.0087798043423074...
     0.0007222414330882...
     0.0000295873218273];
b = [1 ...
     0.0626028287472659...
     0.0238723363198067...
     0.0010727527758408...
     0.0000687107172921];

x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^2/2 + 2*I_precesion_minus(1) - I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % ��� ��������� ������� 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)) + 10^-17);
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
    % ��� ��������� ������� 10^-17
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)-1) + 10^-17);
end

disp('delta = ');
disp((delta)');
grid on, hold on
plot(x, delta, '*')
plot(x_minus, delta_minus, '*')