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

y_star_inv = 1/(sqrt(y_star));

a_last = 1; b_last = 1;
a = zeros(1,N+1); b = zeros(1,N);

link_coefficient = (k+1)*(pi^2)/3;

% задаем лин-триг сетку
baseSize = 2*N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./(y0_inv(end:-1:1)).^2;
x0 = log(exp(y0)-1);

Y = [3.04858735157374161
3.34584190059340303
3.68879069540422755
4.08730270958917252
4.55406258691879984
5.10559265799893147
5.76373546156910699
6.55785012516307209
7.52814427633516026
8.73086555125261121
10.24664082056730230
12.19434940629497000
14.75516278161691552
18.21625034767519935
23.05494184627642795
30.11257710534064103
40.98656328226920920
59.02065112646766210
92.21976738510571181
163.94625312907683679
368.87906954042284724
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797582
28.62575990670503501
]';

z = (I_base*(k+1)./y0).^(2/k);

I_add = [3.97698535404797582
4.50120850382333959
5.13450058508303186
5.90772936928003300
6.86276443443159323
8.05724321262782261
9.57183247081355937
11.52157721162731896
14.07415991914591302
17.48024774155258854
22.12577980514954845
28.62575990670503501
38.00046670694272422
52.02513059266253492
73.97130757221292185
110.31172432162348684
175.06123446242764885
302.39078823846733712
590.48351450010818553
1399.52583034048893751
4723.22570250459648378
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