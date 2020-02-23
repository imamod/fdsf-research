clc 
clear all 
close all
format long

N = 3;
baseSize = 2*N;
x_star = 10;
y_star = log(1+exp(x_star));

for power = 2:4
disp(power);
disp('---------------------');
y_star_inv = 1/(y_star^power);

j = 1:baseSize;
alpha = 2/(2+pi);
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize + (1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./(y0_inv(end:-1:1)).^(1/power);
disp(y0)
end