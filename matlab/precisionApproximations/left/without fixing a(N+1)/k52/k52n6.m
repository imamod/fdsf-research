%% В этом файле значения z и I правильно аппроксимируются, т.е 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 5/2;

N = 6;
x_star = 3;
y_star = log(1+exp(x_star));

a0 = 1; b0 = 1;
a = zeros(1,N+1); b = zeros(1,N);

% задаем лин-триг сетку
baseSize = 2*N+1;
j = 1:baseSize;
alpha = 2/(2+pi);
y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
x0 = log(exp(y0)-1);

f = fopen('D:\repositories\OpenFFD\testOpenFFD\Y_k2.5.txt','r');
Y = fscanf(f, '%f');
fclose(f);

X = log(exp(Y)-1); 

f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_base_k2.5.txt','r');
I_base = fscanf(f, '%f');
fclose(f);

I_base = transpose(I_base);

z = (I_base./(gamma(k+1)*y0)).^(1/k);

f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_add_k2.5.txt','r');
I_add = fscanf(f,'%f');
fclose(f);

I = (I_add./(gamma(k+1)*Y)).^(1/k);

% Задаем матрицы A и B
% B = z(1,:) - (gamma(k+1))*y0.*ones(1,2*N);
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j>=1 && j<=N+1)
            A(i,j) = y0(i)^j;
        elseif (j >= N+2 && j <=baseSize)
            A(i,j) = -z(i)*y0(i)^(j-N-1);
        end
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    if (j>=1 && j<=N+1)
        a(j) = E(j);
    elseif (j >= N+2 && j <=baseSize)
        b(j-N-1) = E(j);
    end
end
% a(N+1) = b(N)*gamma(k+2)^(-1/k);
disp('lg(cond(A)):'); disp(log10(cond(A)));
disp('--------------------------------');
disp('Коэффициенты а:'); disp((10*a)');
disp('----------------------------');
disp('Коэффициенты b:'); disp((10*b)');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*y0(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*y0(j).^m;
    end
    F_base(j) = (a0 + S1)/(b0 + S2);
%    F_base(j) = gamma(k+1)*y0(j)*((a0 + S1)/(b0 + S2))^k;
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% Добавим вспомогательную сетку

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*Y(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*Y(j).^m;
    end
    F(j) = (a0 + S1)/(b0 + S2);
%     F(j) = gamma(k+1)*Y(j)*((a0 + S1)/(b0 + S2))^k;
    delta_additional(j) = F(j)/I(j)-1;
end

disp('lg(dc):');
disp(log10(max(abs(delta_additional))));

grid on, hold on
coeff = 10^15;
xlabel('y'); ylabel('\delta*10^1^5');
plot(Y, delta_additional*coeff, 'k','linewidth', 2.5);
plot(y0, delta_base*coeff, 'k*','linewidth',5)
axis([0 y_star -0.8 0.8]);
line([0;y_star],[0; 0],'linewidth', 2, 'color', 'black');
line([0;log(1+exp(x_star))],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;0],[-1.25; 1.35],'linewidth', 3, 'color', 'black');
% title({'Линейно-тригонометрическая сетка';['k = ', num2str(k), ', N = ', num2str(N), ', x_d_i_v = ', num2str(x_div)]})
% plot(log(1+exp(x_div)),2*10^-16,'b*');