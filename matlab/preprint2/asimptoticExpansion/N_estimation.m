clc
clear all
close all

% ќценка N дл€ различных точностей
format long

for k = -3/2:7/2;
for N = 4:16;
    x = 4*N - 2*k-3;
    p = 1:2*N;
    production = prod(k + 2 - p);
%     disp(production);
    epsilon = production / (x^(2*N));
    output = sprintf('k = %f :N = %d : x = %0.5f : e = %e', k, N, x, epsilon);
    disp(output);
end
    disp('-------------------------------------');
end
