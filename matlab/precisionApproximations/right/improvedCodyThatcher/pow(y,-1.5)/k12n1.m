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

y_star_inv = 1/y_star^(3/2);

a0 = 1; b0 = 1;
a = zeros(1,N); b = zeros(1,N);

% задаем лин-триг сетку
baseSize = 2*N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize +(1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./(y0_inv(end:-1:1)).^(2/3);
x0 = log(exp(y0)-1);

Y = [3.04858735157374161
3.14461562118559801
3.24858150081537733
3.36158964043951602
3.48496772673719102
3.62032714309663417
3.76964490481642489
3.93537622841677415
4.12061216365851113
4.32930512752818508
4.56659958947111733
4.83933076890510616
5.15680169179876291
5.53204143549784177
5.98393828745215650
6.54106408334447842
7.24902499224400021
8.18591243035274552
9.49890993241433712
11.50710989843189402
15.07857961926569779
23.93575314980862245
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797671
7.47184943859957507
]';

 C1 = (k+1)*k*(pi^2)/6;
 z = ((I_base*(k+1)./(y0.^(k+1)) - 1).*(y0.^2))/C1;
%z = (I_base*(k+1)./y0).^(2/k);

I_add = [3.97698535404797582
4.14379591792138058
4.32713428178599013
4.52962504197043536
4.75447999391561549
5.00567170669138672
5.28817238439686488
5.60828852487095286
5.97413932700660588
6.39635631084969347
6.88913325563320456
7.47184943859957684
8.17166749406208126
9.02786382570683088
10.09940720544858905
11.47903293633395272
13.32138927881034363
15.90495580312622792
19.78683840814591832
26.26719265471894360
39.24721570074133581
78.23750666572101409
]';

I = ((I_add*(k+1)./(Y.^(k+1)) - 1).*(Y.^2))/C1;
% grid on, hold on
plot(X, I, 'k*-');

% Задаем матрицы A и B
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j > 0 && j < N + 1)
            A(i,j) = y0(i)^(-2*j);
        elseif (j >= N+1 && j <= baseSize)
            A(i,j) = -z(i)*y0(i)^(-2*(j-N));
        end
    end
end
disp('A matrix: ')
disp(A); 
disp('B matrix: ')
disp(B');
E = A\B';
% disp(E);

for j = 1:length(E)
    if (j > 0 && j <= N)
        a(j) = E(j);
    elseif (j >= N+1 && j <= baseSize)
        b(j-N) = E(j);
    end
end
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
    for n = 1:N
        S1 = S1 + a(n).*y0(j).^(-2*n);
    end
    for m = 1:N
        S2 = S2 + b(m).*y0(j).^(-2*m);
    end
    F_base(j) = (1 + S1)/(1 + S2);
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% Добавим вспомогательную сетку

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0; S2 = 0; 
    for n = 1:N
        S1 = S1 + a(n).*Y(j).^(-2*n);
    end
    for m = 1:N
        S2 = S2 + b(m).*Y(j).^(-2*m);
    end
    F(j) = (1 + S1)/(1 + S2);
    delta_additional(j) = F(j)/I(j)-1;
end

disp('-------------------');
disp('Экстремумы:');
disp(max(abs(delta_additional(1:11))));
disp(max(abs(delta_additional(11:end))));
disp('-------------------');

figure
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