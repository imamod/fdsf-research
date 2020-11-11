%%% Сводные графики по примеру U=exp(x) для аппроксимации многочленом
clc
clear all
close all

% Число узлов
N = 3:12;

delta_start = [
     0.0158306282354284 
 0.00195815749702966
0.000209188606836497
1.95685179664729e-05
1.62667077696632e-06
1.21615331361014e-07
8.25909829416105e-09
5.14106535121073e-10
2.96100921559628e-11
1.58795199212136e-12
];

delta_final = [
    0.00711660918454826
0.000660129126579978
5.25646189293294e-05
3.64183596934797e-06
2.22802543525802e-07
1.21561944954607e-08
6.00109073545241e-10
2.69984035128346e-11
1.11577413974828e-12
4.36317648677687e-14
    ];

log_delta_start = log10(delta_start);
log_delta_final = log10(delta_final);
figure, hold on
xlabel('N'); ylabel('lg(max|P/u|)');
plot(N, log_delta_start, 'k.-', 'MarkerSize', 20);
plot(N, log_delta_final, 'k-o');
axis([N(1) N(end) log_delta_final(end)-0.1 log_delta_start(1)+0.1]);

%% Число итераци
iterations = [6 8 10 12 14 17 19 21 23 22];

figure, hold on
plot(N, iterations, 'k.-', 'MarkerSize', 20);
xlabel('N');
