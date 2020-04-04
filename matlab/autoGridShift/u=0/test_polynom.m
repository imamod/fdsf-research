clc
clear all
close all

format long
%N = 15; % Рабочее N, для которого строим графики 1 - 6
N_arr=15:-3:3;
iterations=[];
L_first_p1 = [];
L_neighb_last = [];
for N = 15:-3:3 % Диагностика зависимости числа итераций от N
% Равномерное распределение точек
n = 0:N;
x = (2.*n-N)/N;

% Максимальное число итераций
maxIterCount = 50;
% Значение константы отношения u_max/u_min - 1
U_RELATION_CONST = 0.01;
L = [];
it = 1;
while it < maxIterCount
    disp('Итерация: '); disp(it);
    % Дополнительные узлы
    clear X; X = [-1];
    for i = 1:length(x)-1
        delta = (x(i+1)-x(i))/20;
        for j = 1:20
            X = [X, x(i) + j*delta];
        end
    end
    % Аппроксимирующий многочлен
    clear U; U = approxPolinom(x, X);

    grid on, hold on
    plot(X, U)
    % Получаем экстремумы
    clear ekstrDelta; ekstrDelta = ekstremList(U, N);
    L = [L, max(abs(ekstrDelta))/min(abs(ekstrDelta)) - 1];
    % Критерий останова
    if L(end) < U_RELATION_CONST
        break;
    end
    x = shiftGridAccurate(x, ekstrDelta); 
    disp('Сдвинутая сетка: '); disp(x');
    disp('-----------------');
    it = it + 1;
end

disp('L = max|e|/min|e| - 1');
disp('L: '); disp(L');
%lg_max_div_min = log10(max_div_min);
%disp('lg(L): '); disp(lg_max_div_min');
disp('-----------------');
neighbours = [];
for i=1:length(L)-1
    neighbours = [ neighbours, L(i)/L(i+1)];
end
disp('Отношение 2 соседних величин L: ');
disp(neighbours');
disp('-----------------');

%%% Сводные диагностики
iterations = [iterations; it];
L_first_p1 = [L_first_p1; L(1) + 1];
L_neighb_last = [L_neighb_last; neighbours(end)];
end

t = table(N_arr', iterations, L_first_p1, L_neighb_last, 'VariableNames',{'N','it', 'L_plus_1', 'L_neighb_last'});
disp(t);