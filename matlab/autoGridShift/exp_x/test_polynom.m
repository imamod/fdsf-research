clc
clear all
close all

format long
%%% Сводные диагностики
%N_arr=3:13; % для аппроксимации многочленом
N_arr=3;%3:10; % для рациональной аппроксимации
iterations=[]; % число итераций
L_first_p1 = []; % диагностика L + 1 на первой итерации
L_neighb_last = []; % диагностика L(n)/L(n+1) на последней итерации
max_delta_first_it = []; % диагностика max|delta| на первой итерации
max_delta_last_it = []; % диагностика max|delta| на последней итерации
% Диагностика зависимости числа итераций от N
for N = N_arr
    % Равномерное распределение точек
    n = 0:N;
    x = (2.*n-N)/N;

    % Максимальное число итераций
    maxIterCount = 50;
    % Значение константы отношения u_max/u_min - 1
    U_RELATION_CONST = 0.02;
    L = [];
    it = 1;
    max_abs_delta_rel = [];
    figure;
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
%         clear U; U = approxPolinom(x, X, exp(x));
        % Рациональная аппроксимация
         clear U; U = rationalApprox(x, X, exp(x));

        % Точное значение экспоненты
        F_prec = exp(X);
        % Относительная погрешность
        delta_rel = U./F_prec - 1;

        grid on, hold on
        plot(X, delta_rel);

        % Получаем экстремумы
        clear ekstrDelta; ekstrDelta = ekstremList(delta_rel, N);
        current_max_abs_delta = max(abs(ekstrDelta));
        if (it > 1 && current_max_abs_delta > max_abs_delta_rel(end))
            disp('break by round error');
            break;
        end
        max_abs_delta_rel = [max_abs_delta_rel, current_max_abs_delta];
        
        disp('max(abs(ekstrDelta)): '); disp(max(abs(ekstrDelta)));
        disp('min(abs(ekstrDelta)): ');disp(min(abs(ekstrDelta)));
        
        L = [L, max(abs(ekstrDelta))/min(abs(ekstrDelta)) - 1];
        disp('L: '); disp(L(end));
        % Критерий останова
        if L(end) < U_RELATION_CONST
            break;
        end
        x = shiftGridAccurate(x, ekstrDelta);
        if isempty(x)
            disp('Ошибка при получении новой сетки');
            break;
        end
        disp('Сдвинутая сетка: '); disp(x');
        disp('-----------------');
        it = it + 1;
    end
    filename = strcat(num2str(N), '.png');
    saveas(gcf, filename);

    disp('max|delta|:'); disp(max_abs_delta_rel');
    disp('-----------------');

    disp('L = max|e|/min|e| - 1');
    disp('L: '); disp(L');
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
    max_delta_first_it = [max_delta_first_it; max_abs_delta_rel(1)];
    max_delta_last_it = [max_delta_last_it; max_abs_delta_rel(end)];
end

%t = table(N_arr', iterations, L_first_p1, L_neighb_last, max_delta_first_it, max_delta_last_it, 'VariableNames',{'N','it', 'L_plus_1', 'L_neighb_last', 'max_delta_1', 'max_delta_end'});
t = table(N_arr', iterations, L_first_p1, max_delta_first_it, max_delta_last_it, 'VariableNames',{'N','it', 'L_plus_1', 'max_delta_1', 'max_delta_end'});
disp(t);