function max_max_delta = commonPart(y0, Y, I_base, I_add, k, N)
baseSize = 2*N;
C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

Z = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% Задаем матрицы A и B
B = (z(1,:) - ones(1,baseSize));
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j > 0 && j < N + 1)
            A(i,j) = (y0(i))^(-2*j);
        elseif (j > N && j < baseSize+1)
            A(i,j) = -z(i)*y0(i)^(2*(N-j));
        end
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

a = zeros(1,N); b = zeros(1,N);
for j = 1:length(E)
    if (j > 0 && j < N + 1)
        a(j) = E(j);
    elseif (j > N && j < baseSize + 1)
        b(j-N) = E(j);
    end
end
disp('lg(cond(A)):'); disp(log10(cond(A)));
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
    F_z_base(j) = (1 + S1)/(1 + S2);
    delta_z(j) = F_z_base(j)/z(j)-1;
%%%     z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);
    F_base(j) = (((F_z_base(j)*2*C1/(k*y0(j)^2)) + 1)^(k/2))*(y0(j)^(k+1))/(k+1);
    delta_base(j) = F_base(j)/I_base(j)-1;
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
    F_z(j) = (1 + S1)/(1 + S2);
    delta_z_additional(j) = F_z(j)/Z(j)-1;
    F(j) = (((F_z(j)*2*C1/(k*Y(j)^2)) + 1)^(k/2))*(Y(j)^(k+1))/(k+1);
    delta_additional(j) = F(j)/I_add(j)-1;
end

figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
title('Относительная погрешность z');
plot(Y,delta_z_additional, 'k','linewidth', 2.5);
plot(y0,delta_z, 'k*','linewidth',5)
% disp((I-F)');

disp('lg(dc):');disp(log10(max(abs(delta_additional))));
disp('-------------------');
disp('Экстремумы:');
span = 11;
max_delta = max(abs(delta_additional(1:span)));
disp('log10(max_delta): '); disp(log10(max_delta));
for i=1:baseSize-1 
    max_delta = max(abs(delta_additional((i*span):(i+1)*span)));
%     disp('Max delta: '); disp(max_delta);
    disp('log10(max_delta): '); disp(log10(max_delta));
end
figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
title('Относительная погрешность I');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)


max_max_delta = max(max_delta);
% disp('max_max_delta: '); disp(max_max_delta');
% disp('-------------------');
end