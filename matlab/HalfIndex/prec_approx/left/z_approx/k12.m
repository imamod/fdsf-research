clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
k = 1/2;

N = 4;
%x_star = 3;
x_star = 0;

y_star = log(1+exp(x_star));

b0 = 1;
a = zeros(1,N-1); b = zeros(1,N);

% задаем лин-триг сетку
baseSize = 2*N-1;
j = 1:baseSize;
alpha = 2/(2+pi);
y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
x0 = log(exp(y0)-1);

%%% global 51
z = [
    0.99947021643942136
0.99751927673517804
0.99558928706854988
0.99481630379222175
0.99475134849058389
    ]';

%%% global 52
% N = 2, x =3
z = [
    0.99771070156454211
1.00182361857576097
1.00336872177306802
    ]';
% x = 0
z = [
    1.00346602055672363
1.00056400789738831
0.99862616353777900
    ]';

% N = 3, x = 0
z =[
    1.00276220477815237
1.00323427708286728
1.00126991590083181
0.99945435829117735
0.99862616353777900
   ]';


% N = 4, x = 0
z = [
    1.00203621225277706
1.00339889116511838
1.00303814782466882
1.00158878581813382
1.00011488278316918
0.99912270810571036
0.99862616353777900
    ]';

% Задаем матрицы A и B
z_tilda= z(1,:) - ones(1,baseSize);
B = z_tilda(1,:);
A = zeros(baseSize,baseSize);
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if (j < N)
            A(i,j) = y0(i)^j;
        elseif (j <= baseSize)
            A(i,j) = -z_tilda(i)*y0(i)^(j-N+1);
        end
    end
end
disp(A);
disp(B');
E = A\B';
disp(E);

for j = 1:length(E)
    if (j < N)
        a(j) = E(j);
    elseif (j <= baseSize)
        b(j-N+1) = E(j);
    end
end

disp('lg(cond(A)):'); disp(log10(cond(A)));
disp('--------------------------------');
disp('Коэффициенты а:'); disp(a');
disp('----------------------------');
disp('Коэффициенты b:'); disp(b');
disp('----------------------------');


disp('Воспроизведение z:');
disp('----------------------------');
z_repr = zeros(1,baseSize);
for j = 1:baseSize
    S1 = 0; S2 = 0; 
    for n = 1:N-1 
        S1 = S1 + a(n).*y0(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*y0(j).^m;
    end
    z_repr(j) = S1/(b0 + S2);
end
disp(z_repr');

disp('Воспроизведение Z:');
disp('----------------------------');
f = fopen('y_add','r');
Y = fscanf(f, '%f'); Y = Y';
fclose(f);

f = fopen('z_add','r');
Z_add = fscanf(f, '%f'); Z_add = Z_add';
fclose(f)

Z_add_tilda = Z_add(1,:) - ones(1,size(Z_add, 2));
Z_repr = zeros(1,size(Z_add_tilda, 2));
for j = 1:size(Z_add_tilda, 2)
    S1 = 0; S2 = 0; 
    for n = 1:N-1
        S1 = S1 + a(n).*Y(j).^n;
    end
    for m = 1:N
        S2 = S2 + b(m).*Y(j).^m;
    end
    Z_repr(j) = S1/(b0 + S2);
end
disp('Погрешность Z:');
disp('----------------------------');
delta_Z = Z_repr - Z_add_tilda;

shift = 20;
for i = 0:baseSize-1
    start = i*shift + 1;
    finish = shift*(i + 1);
    disp(max(abs(delta_Z(start:finish))))
end
