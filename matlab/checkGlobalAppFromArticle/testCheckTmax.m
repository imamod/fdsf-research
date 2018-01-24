clc
clear all
close all

format shortE

k = 1; N = 4;
n = 0:N+1; m = 0:N;

x = [-10 0 10];
y = log(1 + exp(x));

a1 = (gamma(k+2))^(-1/k);
i = 1;
T = x + 39;
while true
    if (i >=4)
        break;
    end
    I_approximate = gamma(k+1)*y.*(1+a1.*y).^k;
    T = x + 39 + log(((T + 1).^k)./I_approximate);
    disp('T for -10    0     10 ')
    disp(T);
    i = i + 1;
end