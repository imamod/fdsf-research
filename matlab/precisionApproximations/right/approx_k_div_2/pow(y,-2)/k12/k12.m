%% z = (Ik/(k!*y))^(1/k)
clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 4;
x_star = 15;
y_star = log(1+exp(x_star));

a_last = 1; b_last = 1;
a = zeros(1,N+1); b = zeros(1,N);
link_coefficient = (k+1)*(pi^2)/3;

% задаем лин-триг сетку
baseSize = 2*N;
j = 1:baseSize;
% y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
f = fopen('n4/y_base_right','r');
y0 = fscanf(f, '%f');
y0 = y0';
fclose(f);
x0 = log(exp(y0)-1);

f = fopen('n4/y_add_right','r');
Y = fscanf(f, '%f');
Y = Y';
fclose(f); 
X = log(exp(Y)-1); 

f = fopen('n4/I_base_right','r');
I_base = fscanf(f, '%f');
I_base = I_base';
fclose(f); 
z = (I_base*(k+1)./y0).^(2/k);

f = fopen('n4/I_add_right','r');
I_add = fscanf(f, '%f');
I_add = I_add';
fclose(f); 
I = (I_add*(k+1)./Y).^(2/k);

% Строим z
grid on, hold on
plot(X, I, 'k*');

% Задаем матрицы A и B
B = (z(1,:) - ones(1,baseSize).*(y0.^2)).*y0.^(2*N) + link_coefficient*y0.^(2*N);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j > 0 && j < N + 1)
            A(i,j) = (y0(i))^(2*(j-1));
        elseif (j >= N+1 && j <=baseSize-1)
            A(i,j) = -z(i)*y0(i)^(2*(j-N-1));
        elseif (j == baseSize)
            A(i,j) = -z(i)*y0(i)^(2*N-2) + (y0(i)^(2*N));
        end
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    if (j > 0 && j <= N)
        a(j) = E(j);
    elseif (j >= N+1 && j <= baseSize)
        b(j-N) = E(j);
    end
end
a(N+1) = b(N) - link_coefficient;
% disp('lg(cond(A)):'); disp(log10(cond(A)));
disp('--------------------------------');
disp('Коэффициенты а:'); disp(a');
disp('----------------------------');
disp('Коэффициенты b:'); disp(b');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*y0(j).^(2*(n-1));
    end
    for m = 1:N
        S2 = S2 + b(m).*y0(j).^(2*(m-1));
    end
    F_base(j) = (y0(j).^(2*N+2) + S1)/(y0(j).^(2*N) + S2);
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% Добавим вспомогательную сетку

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*Y(j).^(2*(n-1));
    end
    for m = 1:N
        S2 = S2 + b(m).*Y(j).^(2*(m-1));
    end
    F(j) = (Y(j).^(2*N+2) + S1)/(Y(j).^(2*N) + S2);
    delta_additional(j) = F(j)/I(j)-1;
end

disp('lg(dc):');disp(log10(max(abs(delta_additional))));

disp('-------------------');
disp('Экстремумы Z:');
span = 11;
max_delta = max(abs(delta_additional(1:span)));
disp(log10(max_delta));
for i=1:baseSize-1 
    max_delta = max(abs(delta_additional((i*span):(i+1)*span)));
    disp(log10(max_delta));
end
disp('-------------------');

figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)

app_I_base = zeros(1,length(F_base));
delta_base_I = zeros(1,length(F_base));
for i=1:length(F_base)
    app_I_base(i) = F_base(i)^(k/2)*(y0(i)/(k+1));
    delta_base_I(i) = app_I_base(i)/I_base(i)-1;
end

app_I = zeros(1,length(Y));
delta_additional_I = zeros(1,length(Y));
for i=1:length(app_I)
    app_I(i) = F(i)^(k/2)*(Y(i)/(k+1));
    delta_additional_I(i) = app_I(i)/I_add(i)-1;
end

disp('lg(dc):');disp(log10(max(abs(delta_additional_I))));
disp('-------------------');
disp('Экстремумы I:');
span = 11;
max_delta = max(abs(delta_additional_I(1:span)));
disp(log10(max_delta));
for i=1:baseSize-1 
    max_delta = max(abs(delta_additional_I((i*span):(i+1)*span)));
    disp(log10(max_delta));
end
disp('-------------------');
