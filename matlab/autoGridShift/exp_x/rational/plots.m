%%% Сводные графики по примеру U=exp(x)
clc
clear all
close all

% Число узлов
N = 3:10;

delta_start = [
     0.00339856813527284 
0.000159364982344368
1.12876570379505e-05
4.39628657700197e-07
2.29408589902391e-08
7.49799777821636e-10
3.06407121897223e-11
8.57425241918008e-13
];

delta_final = [
    0.00235685332445557
0.000110933041667094
5.33294895876946e-06
1.85058270263383e-07
6.47220721461395e-09
1.76745396096578e-10
4.84812190393313e-12
1.09356967925578e-13
    ];

log_delta_start = log10(delta_start);
log_delta_final = log10(delta_final);
figure, hold on
xlabel('n'); ylabel('lg(|p_n_-_1_/_2|)');
plot(N, log_delta_start, 'k.-', 'MarkerSize', 20);
plot(N, log_delta_final, 'k-o');
axis([N(1) N(end) log_delta_final(end)-0.1 log_delta_start(1)+0.1]);

%% Число итераци
iterations = [5 4 8 6 10 8 12 11];

figure, hold on
plot(N, iterations, 'k.-', 'MarkerSize', 20);
xlabel('n'); ylabel('it');

%% СРавнение полинома и многочлена

relation_pol_rat = [3.02 5.95 9.85 19.67 34.43 68.78 123.78 246.88];
figure, hold on
plot(3:10, log10(relation_pol_rat), 'k.-', 'MarkerSize', 20);
xlabel('n'); ylabel('lg(\delta_p/\delta_r)');
