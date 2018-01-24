clc
clear all
close all

k = 1/2;

N = [3, 5, 7, 9, 11, 13];

max_max_delta = [];
for j = 1:length(N)
%     delta_additional = readData(k , N(j));
%     span = 11;
%     max_delta = max(abs(delta_additional(1:span)));
%     for i=1:2*N-1 
%         max_delta = [max_delta, max(abs(delta_additional((i*span):(i+1)*span)))];
%     end
%     disp('max_delta: '); disp(max_delta');
%     max_max_delta = [max_max_delta, max(max_delta)];
	max_max_delta = [max_max_delta, readData(k , N(j))];
end
disp('max_max_delta: ');
disp(max_max_delta')

figure
grid on, hold on
plot(N, log10(max_max_delta),'k*-', 'linewidth', 3)
xlabel('N'); ylabel('log_1_0\delta');

