clc
clear all
close all

format long
%%% ������� �����������
%N_arr=3:13; % ��� ������������� �����������
N_arr=3;%3:10; % ��� ������������ �������������
iterations=[]; % ����� ��������
L_first_p1 = []; % ����������� L + 1 �� ������ ��������
L_neighb_last = []; % ����������� L(n)/L(n+1) �� ��������� ��������
max_delta_first_it = []; % ����������� max|delta| �� ������ ��������
max_delta_last_it = []; % ����������� max|delta| �� ��������� ��������
% ����������� ����������� ����� �������� �� N
for N = N_arr
    % ����������� ������������� �����
    n = 0:N;
    x = (2.*n-N)/N;

    % ������������ ����� ��������
    maxIterCount = 50;
    % �������� ��������� ��������� u_max/u_min - 1
    U_RELATION_CONST = 0.02;
    L = [];
    it = 1;
    max_abs_delta_rel = [];
    figure;
    while it < maxIterCount
        disp('��������: '); disp(it);
        % �������������� ����
        clear X; X = [-1];
        for i = 1:length(x)-1
            delta = (x(i+1)-x(i))/20;
            for j = 1:20
                X = [X, x(i) + j*delta];
            end
        end

        % ���������������� ���������
%         clear U; U = approxPolinom(x, X, exp(x));
        % ������������ �������������
         clear U; U = rationalApprox(x, X, exp(x));

        % ������ �������� ����������
        F_prec = exp(X);
        % ������������� �����������
        delta_rel = U./F_prec - 1;

        grid on, hold on
        plot(X, delta_rel);

        % �������� ����������
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
        % �������� ��������
        if L(end) < U_RELATION_CONST
            break;
        end
        x = shiftGridAccurate(x, ekstrDelta);
        if isempty(x)
            disp('������ ��� ��������� ����� �����');
            break;
        end
        disp('��������� �����: '); disp(x');
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
    disp('��������� 2 �������� ������� L: ');
    disp(neighbours');
    disp('-----------------');

    %%% ������� �����������
    iterations = [iterations; it];
    L_first_p1 = [L_first_p1; L(1) + 1];
    L_neighb_last = [L_neighb_last; neighbours(end)];
    max_delta_first_it = [max_delta_first_it; max_abs_delta_rel(1)];
    max_delta_last_it = [max_delta_last_it; max_abs_delta_rel(end)];
end

%t = table(N_arr', iterations, L_first_p1, L_neighb_last, max_delta_first_it, max_delta_last_it, 'VariableNames',{'N','it', 'L_plus_1', 'L_neighb_last', 'max_delta_1', 'max_delta_end'});
t = table(N_arr', iterations, L_first_p1, max_delta_first_it, max_delta_last_it, 'VariableNames',{'N','it', 'L_plus_1', 'max_delta_1', 'max_delta_end'});
disp(t);