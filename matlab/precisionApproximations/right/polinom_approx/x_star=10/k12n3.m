%% � ���� ����� �������� z � I ��������� ����������������, �.� 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% ����������� k ������� �����-������
k = 1/2;

N = 3;

x_star = 10;
y_star = log(1+exp(x_star));

y_star_inv = 1/y_star^2;

% ������ ���-���� �����
baseSize = N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize +(1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);

Y = [10.00004539889921951
10.13094257365187545
10.26711828953063765
10.40893712872289179
10.55679993000528150
10.71114856026497897
10.87247147667029168
11.04131023913890175
11.21826717114195837
11.40401441608143429
11.59930470000498026
11.80498419410607447
12.14449633116664451
12.51509424418397032
12.92183089635294557
13.37098935895881624
13.87049671072096402
14.43052195699585738
15.06436319666928192
15.78980542024332223
16.63127564590926610
17.62341524776945789
18.81732134693724490
19.73577312752885859
20.80333148911197583
22.06526515083407602
23.58876066935541260
25.47877354914904302
27.91059802086978792
31.20499723366796729
36.03242710583946007
44.13053030166815205
62.40999446733593459
]';

X = log(exp(Y)-1); 

I_base = [21.34447149235519703
27.28097188170979948
54.60843505826677813
]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

I_add = [21.34447149235520413
21.75799837907527845
22.19109445995611551
22.64525263686044099
23.12212547853830458
23.62354748380591118
24.15156121605255279
24.70844812396543944
25.29676506664164748
25.91938782369225081
26.57956321237986330
27.28097188170978527
28.45234732536453492
29.75001520076441253
31.19674648040383502
32.82130251007851740
34.66058179022784458
36.76277415521101233
39.19212668793978338
42.03638264047911122
45.41884388961057084
49.51884023572544180
54.60843505826677813
58.63628938148695369
63.43761606335392145
69.27451984058841106
76.54728361275184056
85.90178646426244313
98.45787769568291026
116.35775656941413558
144.33171437415717264
195.56566615869840575
328.79675585154780038
]';

I = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% ������ ������� A � B
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        A(i,j) = y0(i)^(-2*j);
    end
end
disp('A matrix: ')
disp(A); 
disp('B matrix: ')
disp(B');
E = A\B';

a0 = 1;
a = zeros(1,baseSize);
for j = 1:length(E)
    a(j) = E(j);
end
disp('--------------------------------');
disp('������������ �:'); disp(a');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    F_base(j) = 1; 
    for n = 1:baseSize
        F_base(j) = F_base(j) + a(n).*y0(j).^(-2*n);
    end
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% ������� ��������������� �����

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    F(j) = 1;
    for n = 1:baseSize
        F(j) = F(j) + a(n).*Y(j).^(-2*n);
    end
    delta_additional(j) = F(j)/I(j)-1;
end

disp('-------------------');
disp('����������:');
disp(max(abs(delta_additional(1:11))));
disp(max(abs(delta_additional(11:22))));
disp(max(abs(delta_additional(22:end))));
disp('-------------------');

grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)
% axis([0 y_star -1.5*10^(-1) 1.5*10^(-1)])
% line([0;y_star],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;log(1+exp(x_star))],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;0],[-1.25; 1.35],'linewidth', 3, 'color', 'black');
% title({'�������-������������������ �����';['k = ', num2str(k), ', N = ', num2str(N), ', x_d_i_v = ', num2str(x_div)]})
% plot(log(1+exp(x_div)),2*10^-16,'b*');