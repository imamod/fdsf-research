clc
clear all
close all

format long
%N = 15; % Рабочее N, для которого строим графики 1 - 6
N_arr=15%3:15;
iterations=[];
L_first_p1 = [];
L_neighb_last = [];
%%%
min_ekstr_it = [];
max_ekstr_it = [];
final_ekstr_it = [];
for N = N_arr % Диагностика зависимости числа итераций от N
% Равномерное распределение точек
n = 0:N;
x = (2.*n-N)/N;

% Максимальное число итераций
maxIterCount = 50;
% Значение константы отношения u_max/u_min - 1
U_RELATION_CONST = 0.01;
L = [];
it = 1;
% Профиль погрешности за 
first_it_profile = [];
first_it_ekstr = [];
first_it_X = [];
last_it_profile = [];
last_it_X = [];
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
    if it == 1
      %  min_ekstr_it = [min_ekstr_it, min(abs(ekstrDelta))];
      %  max_ekstr_it = [max_ekstr_it, max(abs(ekstrDelta))];
        first_it_profile = U;
        first_it_ekstr = ekstrDelta;
        first_it_X = X;
    end
    % Критерий останова
    if L(end) < U_RELATION_CONST
        last_it_profile = U;
        last_it_X = X;
      %  final_ekstr_it = [final_ekstr_it, max(abs(ekstrDelta))];
        break;
    end
    x = shiftGridAccurate(x, ekstrDelta); 
    disp('Сдвинутая сетка: '); disp(x');
    disp('-----------------');
    it = it + 1;
end

disp('L = max|e|/min|e| - 1');
disp('L: '); disp(L');
lg_L = log10(L);
disp('lg(L): '); disp(lg_L');
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

if not isempty(min_ekstr_it)
    ekstr = table(N_arr', min_ekstr_it', max_ekstr_it', final_ekstr_it', 'VariableNames',{'N', 'min_ekstr_it', 'max_ekstr_it', 'final_ekstr_it'});
    disp(ekstr);
end