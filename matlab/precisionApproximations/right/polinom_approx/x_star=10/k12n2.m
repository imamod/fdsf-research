%% В этом файле значения z и I правильно аппроксимируются, т.е 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 2;

x_star = 10;
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

Y = [10.00004539889921773
10.23537278163795960
10.48813609646871114
10.76060058869978953
11.05546615822286682
11.37598082535665789
11.72609263448666717
12.11065639718728626
12.53572032109215684
13.00893177060993722
13.54012547812673262
14.14219982746994120
14.83246431164086587
15.63479017931555504
16.58319923733430201
17.72818569220453710
19.14862908740031244
20.97627219293742584
23.45218526897333433
27.08025095625346523
33.16639847466860402
46.90437053794666866
]';

X = log(exp(Y)-1); 

I_base = [21.34447149235519703
35.67524537583092581
]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

I_add = [21.34447149235519703
22.08986770916970954
22.90025254993765813
23.78497910058832687
24.75533660646742717
25.82508779014026246
27.01119788172135472
28.33484097643458099
29.82281634560274952
31.50958570949543258
33.44027748846910697
35.67524537583092581
38.29721876297616490
41.42296517503770303
45.22322004296682252
49.95873397637950575
56.05025216060314364
64.22729470033510779
75.88541803308312694
94.10628717497553453
127.48039525604447419
214.27538056783359366
]';

I = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% Задаем матрицы A и B
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
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
a = zeros(1,baseSize);
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
    for n = 1:baseSize
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
disp(max(abs(delta_additional(11:22))));
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