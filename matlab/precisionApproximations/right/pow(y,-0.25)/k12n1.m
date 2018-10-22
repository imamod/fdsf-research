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
% x_star = 2;
x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/((y_star)^(1/4));

a_last = 1; b_last = 1;
a = zeros(1,N+1); b = zeros(1,N);

link_coefficient = (k+1)*(pi^2)/3;

% задаем лин-триг сетку
baseSize = 2*N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./(y0_inv(end:-1:1)).^4;
x0 = log(exp(y0)-1);

Y = [3.04858735157374738
3.67208045325898036
4.46343674143912317
5.47992939457385742
6.80298238292809376
8.55054272135462234
10.89706235702910497
14.10666426923969041
18.58990729462357550
25.00437234796607555
34.44009831357345064
48.77739762517992261
71.41498786302591384
108.84771812684941494
174.35299771246550904
297.43851671397692371
551.04157301717486916
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797582
227.22802709098678520
]';

z = (I_base*(k+1)./y0).^(2/k);

I_add = [3.97698535404798559
5.10294699540127805
6.67325654367720311
8.90711093383547947
12.14923759362519462
16.95324804856643652
24.23239905671557182
35.54197571239791387
53.62609486755162180
83.51990569491421468
134.88297210479089472
227.22802709098698415
402.43695176532628466
757.15179786534019968
1534.86610129024370508
3419.87795076215843437
8623.56921035348750593
]';


I = (I_add*(k+1)./Y).^(2/k);

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

disp('-------------------');
disp('Экстремумы:');
disp(max(abs(delta_additional(1:11))));
disp(max(abs(delta_additional(11:end))));
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