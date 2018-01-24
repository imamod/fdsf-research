clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 1;
baseSize = 2*N;
% x_star = 2;
x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/(y_star^2);

% a0 = 1; b0 = 1;
a = zeros(1,N); b = zeros(1,N);

% задаем лин-триг сетку 
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);

Y = [3.04858735157374205
3.12032863412584449
3.19738538875014910
3.28044819213102778
3.37034012859779519
3.46805138102125898
3.57478553975496238
3.69202262979574991
3.82160649170195033
3.96586848073037856
4.12780678734541873
4.31135357867466062
4.52178578090403160
4.76638071967328347
5.05551019289669323
5.40456773061796181
5.83760034151960561
6.39477077750029821
7.14957107950992654
8.25561357469083745
10.11102038579338647
14.29914215901985308
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797671
6.35962395754328824
]';

% z = (I_base*(k+1)./y0).^(2/k);
C1 = (k*pi^2)/12;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)/C1;

I_add = [3.97698535404797671
4.10137679619248274
4.23649792045401519
4.38389672313587475
4.54544297334555480
4.72341783408510185
4.92063554805746417
5.14061153583199371
5.38779914291691764
5.66793042998608954
5.98851907086318480
6.35962395754328735
6.79504774608824658
7.31429289559311968
7.94590618981402130
8.73353104508873024
9.74766204086129306
11.11064459779278302
13.05672278606514247
16.10352981695412566
21.69488436865299263
36.26584742434106090
]';

% I = (I_add*(k+1)./Y).^(2/k);
I = ((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2)/C1;

% Задаем матрицы A и B
B = (z(1,:) - ones(1,baseSize));
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j > 0 && j < N + 1)
            A(i,j) = (y0(i))^(-2*j);
        elseif (j > N && j < baseSize+1)
            A(i,j) = -z(i)*y0(i)^(2*(N-j));
        end
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    if (j > 0 && j < N + 1)
        a(j) = E(j);
    elseif (j > N && j < baseSize + 1)
        b(j-N) = E(j);
    end
end
disp('lg(cond(A)):'); disp(log10(cond(A)));
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

disp((I-F)');

disp('lg(dc):');disp(log10(max(abs(delta_additional))));
disp('-------------------');
disp('Экстремумы:');
disp(max(abs(delta_additional(1:11))));
disp(max(abs(delta_additional(11:22))));
disp('-------------------');

figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)
