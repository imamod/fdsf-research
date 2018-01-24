clc
clear all
close all

k = 1/2;

x_star = 3;
% x_star = 5;
y_star = log(1+exp(x_star));
y_star_inv = 1/(y_star);

% N = [3, 5, 7, 9]; %, 10, 11, 13];
N = [9];

max_max_delta = [];
for j = 1:length(N)
    baseSize = 2*N;

    % задаем лин-триг сетку 
    j = 1:baseSize;
    alpha = 2/(2+pi);
    % alpha = 0.7;
    % y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
    % y0_inv = 0.5*y_star_inv*(2*(1-alpha)*j/(baseSize)+alpha*(1 - cos(pi*j/(baseSize))));
    % y0_inv = j * y_star_inv / baseSize;
    % y0 = 1./(y0_inv(end:-1:1));
    % x0 = log(exp(y0)-1);

    %% Read before shift calculated values
%     prefix = strcat('testShift/tau075/nonshifted/k12/n', num2str(N));
  
    % prefix = strcat('k12/n', num2str(N));
    % prefix = strcat('linearGrid/k12/n', num2str(N));
    % prefix = strcat('linTrigAlpha07/k12/n', num2str(N));
    % prefix = strcat('changedWeights/k12/n', num2str(N));
    % prefix = strcat('testShift/k12/n', num2str(N));
    %% Read after shift recalculated values
    prefix = strcat('testShift/tau1/k12/n', num2str(N));
    % prefix = strcat('x5linTrig/k12/n', num2str(N));

    filename = strcat(prefix, '/y0.txt');
    disp(filename)
    y0 = readFile(filename);

    filename = strcat(prefix, '/Y.txt');
    Y = readFile(filename);

    filename = strcat(prefix, '/I_base.txt');
    I_base = readFile(filename);

    filename = strcat(prefix, '/I_add.txt');
    I_add = readFile(filename);

    max_max_delta = [max_max_delta, commonPart(y0, Y, I_base, I_add, k, N)];
    % 	max_max_delta = [max_max_delta, readData(k , N(j))];
    
end
% disp('max_max_delta: ');
% disp(max_max_delta')
% 
% figure
% grid on, hold on
% plot(N, log10(max_max_delta),'k*-', 'linewidth', 3)
% xlabel('N'); ylabel('log_1_0\delta');

