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

y_star_inv = 1/y_star;

a0 = 1; b0 = 1;
a = zeros(1,N); b = zeros(1,N);

% задаем лин-триг сетку
baseSize = 2*N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize +(1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./y0_inv(end:-1:1);
x0 = log(exp(y0)-1);

Y = [3.04858735157374205
3.19375817783915839
3.35344608673111599
3.52994324919064839
3.72605120747901752
3.94523069027190054
4.19180760841389422
4.47126144897482103
4.79063726675873713
5.15914782574017927
5.58907681121852740
6.09717470314748500
6.70689217346223288
7.45210241495803594
8.38361521682779021
9.58127453351747604
11.17815362243705479
13.41378434692446575
16.76723043365558397
22.35630724487411314
33.53446086731116793
67.06892173462233586
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797671
10.37466055155188727
]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1)) - 1).*(y0.^2))/C1;

I_add = [3.97698535404797671
4.23010247261796124
4.51492209355336094
4.83742261873394419
5.20511314208465858
5.62753834617653226
6.11699164947886853
6.68954211977156721
7.36654562606189245
8.17692390105180955
9.16069965132320085
10.37466055155188727
11.90178355856969183
13.86763084224593179
16.47043781672759621
20.04003589567702548
25.16298675926659101
32.97754481359998380
45.97358609901864668
70.64499610295428056
129.60515765677223499
366.27716559854576417
]';

I = ((I_add*(k+1)./(Y.^(k+1)) - 1).*(Y.^2))/C1;
grid on, hold on
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