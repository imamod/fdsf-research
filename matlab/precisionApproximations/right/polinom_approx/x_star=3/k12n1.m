%% В этом файле значения z и I правильно аппроксимируются, т.е 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 1;

x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/y_star^2;

% задаем лин-триг сетку
baseSize = N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize +(1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);

Y = [3.04858735157374205
3.19738538875014910
3.37034012859779564
3.57478553975496327
3.82160649170195077
4.12780678734541873
4.52178578090403160
5.05551019289669323
5.83760034151960561
7.14957107950992654
10.11102038579338647
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797671]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

I_add = [3.97698535404797671
4.23649792045401519
4.54544297334555480
4.92063554805746595
5.38779914291691941
5.98851907086318480
6.79504774608824658
7.94590618981402130
9.74766204086129306
13.05672278606514247
21.69488436865299263
]';

I = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% Задаем матрицы A и B
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:baseSize
    for j = 1:baseSize
        A(i,j) = y0(i)^(-2*j);
    end
end
disp('A matrix: ')
disp(A); 
disp('B matrix: ')
disp(B');
E = A\B';
% disp(E);

a0 = 1;
a = zeros(1,N);
for j = 1:length(E)
    a(j) = E(j);
end
disp('--------------------------------');
disp('Коэффициенты а:'); disp(a');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    S1 = 0; 
    for n = 1:N
        S1 = S1 + a(n).*y0(j).^(-2*n);
    end
    F_base(j) = 1 + S1;
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% Добавим вспомогательную сетку

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0;
    for n = 1:baseSize
        S1 = S1 + a(n).*Y(j).^(-2*n);
    end
    F(j) = 1 + S1;
    delta_additional(j) = F(j)/I(j)-1;
end

disp('-------------------');
disp('Экстремумы:');
disp(max(abs(delta_additional(1:11))));
disp('-------------------');

grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)
% axis([0 y_star -1.5*10^(-1) 1.5*10^(-1)])
% line([0;y_star],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;log(1+exp(x_star))],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;0],[-1.25; 1.35],'linewidth', 3, 'color', 'black');
% title({'Линейно-тригонометрическая сетка';['k = ', num2str(k), ', N = ', num2str(N), ', x_d_i_v = ', num2str(x_div)]})
% plot(log(1+exp(x_div)),2*10^-16,'b*');