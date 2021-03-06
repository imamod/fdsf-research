%% � ���� ����� �������� z � I ��������� ����������������, �.� 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% ����������� k ������� �����-������
k = 3/2;

N = 2;
% x_star = 2;
x_star = 3;
y_star = log(1+exp(x_star));

a0 = 1; b0 = 1;
a = zeros(1,N+1); b = zeros(1,N);

% ������ ���-���� ����� 
j = 1:2*N;
alpha = 2/(2+pi);
y0 = 0.5*y_star*(2*alpha*j/(2*N)+(1-alpha)*(1 - cos(pi*j/(2*N))));
x0 = log(exp(y0)-1);

Y = [0.05175040116767015
0.10350080233534030
0.15525120350301044
0.20700160467068060
0.25875200583835073
0.31050240700602089
0.36225280817369104
0.41400320934136120
0.46575361050903136
0.51750401167670146
0.56925441284437162
0.65607616402096247
0.74289791519755333
0.82971966637414418
0.91654141755073504
1.00336316872732589
1.09018491990391664
1.17700667108050738
1.26382842225709813
1.35065017343368887
1.43747192461027962
1.52429367578687036
1.61111542696346133
1.69793717814005229
1.78475892931664326
1.87158068049323423
1.95840243166982519
2.04522418284641594
2.13204593402300668
2.21886768519959743
2.30568943637618817
2.39251118755277892
2.47933293872936966
2.53108333989703960
2.58283374106470953
2.63458414223237947
2.68633454340004940
2.73808494456771934
2.78983534573538927
2.84158574690305921
2.89333614807072914
2.94508654923839908
2.99683695040606901
3.04858735157373895
]';

X = log(exp(Y)-1); 

I_base = [0.90968755107532462
3.30055373767231464
7.16291360777494202
10.35371486476145009
]';

z = (I_base./(gamma(k+1)*y0)).^(1/k);

I_add = [0.06995459910770487
0.14227082334439811
0.21700959692607633
0.29423261438840997
0.37400232160856928
0.45638189603384216
0.54143522615253337
0.62922689024469081
0.71982213445212084
0.81328685020891578
0.90968755107532462
1.07819905696548179
1.25548406991198691
1.44186730136312846
1.63767734252099650
1.84324631942080797
2.05890954319919128
2.28500515738879706
2.52187378405642360
2.76985817056394623
3.02930283867652506
3.30055373767231419
3.58395790302389106
3.87986312212512763
4.18861760843047470
4.51056968525823176
4.84606748038730384
5.19545863244997719
5.55909000999369152
5.93730744395348164
6.33045547414661947
6.73887711027236769
7.16291360777493491
7.42322922321632017
7.68928455843957614
7.96115081943457259
8.23889891305797306
8.52259943472901682
8.81232265679960491
9.10813851759262327
9.41011661110102793
9.71832617733906545
10.03283609333599102
10.35371486476143232
]';

I = (I_add./(gamma(k+1)*Y)).^(1/k);

% ������ ������� A � B
B = z(1,:) - ones(1,2*N);
A = zeros(2*N,2*N);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j>=1 && j<=N)
            A(i,j) = (y0(i))^(j);
        elseif (j >= N+1 && j <=2*N-1)
            A(i,j) = -z(i)*y0(i)^(j-N);
        elseif (j == 2*N)
            A(i,j) = -z(i)*y0(i)^N + (gamma(k+2)^(-1/k))*(y0(i)^(N+1));
        end
    end
end
disp(A); 
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    if (j>=1 && j<=N)
        a(j) = E(j);
    elseif (j >= N+1 && j <=2*N)
        b(j-N) = E(j);
    end
end
a(N+1) = b(N)*gamma(k+2)^(-1/k);
disp('lg(cond(A)):'); disp(log10(cond(A)));
disp('--------------------------------');
disp('������������ �:'); disp(a);
disp('----------------------------');
disp('������������ b:'); disp(b);
disp('----------------------------');

F_base = zeros(1,2*N);
delta_base = zeros(1,2*N);
for j = 1:2*N
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*y0(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*y0(j).^m;
    end
    F_base(j) = (a0 + S1)/(b0 + S2);
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% ������� ��������������� �����

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    S1 = 0; S2 = 0; 
    for n = 1:N+1 
        S1 = S1 + a(n).*Y(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*Y(j).^m;
    end
    F(j) = (a0 + S1)/(b0 + S2);
    delta_additional(j) = F(j)/I(j)-1;
end

disp('lg(dc):');disp(log10(max(abs(delta_additional))));

grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)
% axis([0 y_star -1.5*10^(-1) 1.5*10^(-1)])
line([0;y_star],[0; 0],'linewidth', 2, 'color', 'black');
line([0;log(1+exp(x_star))],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;0],[-1.25; 1.35],'linewidth', 3, 'color', 'black');
% title({'�������-������������������ �����';['k = ', num2str(k), ', N = ', num2str(N), ', x_d_i_v = ', num2str(x_div)]})
% plot(log(1+exp(x_div)),2*10^-16,'b*');