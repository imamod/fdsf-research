clc
clear all
close all

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

x = -5:0.01:20;
y = log(1 + exp(x));
for i = 1:length(y)
    I_precesion(i) = gamma(k+1)*y(i).*(sum(a.*(y(i).^n))/sum(b.*(y(i).^m))).^k;

    I_approximate1(i) = y(i).*((gamma(k+1))^(1/k) + y(i)/((k+1)^(1/k))).^k;
    I_approximate2(i) = y(i).*sqrt(1 + y(i)*y(i)/4);
    delta1(i) = log(I_approximate1(i)/I_precesion(i));
    delta2(i) = log(I_approximate2(i)/I_precesion(i));
end

grid on, hold on
plot(x, delta1, '*')
plot(x, delta2, 'k*')