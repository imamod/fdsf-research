clc 
clear all 
close all
format long

% ����������� k ������� �����-������
k = 1/2;

N = 1;
baseSize = N;
x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/(y_star^2);

a = zeros(1,N); b = zeros(1,N);

% ������ ���-���� ����� 
j = 1:baseSize;
y0_inv = y_star_inv*((sin(pi*j/(2*baseSize))).^2);
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);

Y = [3.04858735157374205
3.19738538875014910
3.37034012859779564
3.57478553975496327
3.82160649170195077
4.12780678734541873
4.52178578090403160
5.05551019289669323
5.83760034151960561
7.14957107950992654
10.11102038579338647
]';

X = log(exp(Y)-1); 

I_base = [3.97698535404797671
]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

I_add = [3.97698535404797671
4.23649792045401519
4.54544297334555480
4.92063554805746595
5.38779914291691941
5.98851907086318480
6.79504774608824658
7.94590618981402130
9.74766204086129306
13.05672278606514247
21.69488436865299263
]';

I = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% ������ ������� A � B
B = (z(1,:) - ones(1,baseSize));
A = zeros(baseSize,baseSize);
for i = 1:baseSize
    for j = 1:baseSize
        A(i,j) = y0(i)^(-2*j);
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    a(j) = E(j);
end
disp('lg(cond(A)):'); disp(log10(cond(A)));
disp('--------------------------------');
disp('������������ �:'); disp(a');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    S1 = 0; 
    for n = 1:N
        S1 = S1 + a(n).*y0(j).^(-2*n);
    end
    F_base(j) = 1 + S1;
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% ������� ��������������� �����

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0;
    for n = 1:baseSize
        S1 = S1 + a(n).*Y(j).^(-2*n);
    end
    F(j) = 1 + S1;
    delta_additional(j) = F(j)/I(j)-1;
end

% disp((I-F)');

disp('lg(dc):');disp(log10(max(abs(delta_additional))));
disp('-------------------');
disp('����������:');
disp(max(abs(delta_additional(1:11))));
disp('-------------------');

figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)

figure
grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,I, 'k','linewidth', 2.5);
plot(y0,z, 'k*','linewidth',5)
