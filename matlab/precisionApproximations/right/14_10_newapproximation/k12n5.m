clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 5;
baseSize = 2*N;
% x_star = 2;
x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/(y_star^2);

a = zeros(1,N); b = zeros(1,N);

% задаем лин-триг сетку 
j = 1:baseSize;
alpha = 2/(2+pi);
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./sqrt(y0_inv(end:-1:1));
% x0 = log(exp(y0)-1);

f = fopen('k12/n5/Y_12.txt','r');
Y = fscanf(f, '%f');
Y = Y';
fclose(f); 

f = fopen('k12/n5/I_base_12.txt','r');
I_base = fscanf(f, '%f');
I_base = I_base';
fclose(f);

f = fopen('k12/n5/I_add_12.txt','r');
I_add = fscanf(f, '%f');
I_add = I_add';
fclose(f);

commonPart(y0, Y, I_base, I_add, k, N);
