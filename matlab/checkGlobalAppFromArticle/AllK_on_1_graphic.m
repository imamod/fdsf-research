%% Проверка глобальной аппроксимации из статьи
clc
clear all
close all

format long

k = 1; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.2715113138214362780964488 ...
     0.0562661238060587633169245 ...
     0.0067420740469345689743873 ...
     0.0005169505155333205859985 ...
     0.0000194771836765773190602];
b = [1 ...
     0.0215113138214352840651584 ...
     0.0231105175729721417901084 ...
     0.0003669081577365413477999 ...
     0.0000610424408732720110769];

x = 0:0.1:30;
y = log(1 + exp(x));

x_minus = 0:-0.1:-30;
y_minus = log(1 + exp(x_minus));

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^2/2 + pi^2/6 - I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % Для страховки добавим 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)));
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
    % Для страховки добавим 10^-17
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)) + 10^-17);
end

grid on, hold on
plot(x, delta, 'k', 'linewidth', 2)
plot(x_minus, delta_minus, 'k', 'linewidth', 2)

%% Проверка глобальной аппроксимации из статьи

k = 2; N = 4;
n = 0:N+1; m = 0:N;

a = [1 ...
     0.2263816364340698560028783 ...
     0.0533684335574798857246766 ...
     0.0062904756340795211604491 ...
     0.0005023228274452983506998 ...
     0.0000189379675088061004880];
 
b = [1 ...
     0.0388816364340691133155655 ...
     0.0243043998742774445085992 ...
     0.0006290985326433190105734 ...
     0.0000657018161945458806177];

x = 0:0.1:30;
y = log(1 + exp(x));

x_minus = 0:-0.1:-30;
y_minus = log(1 + exp(x_minus));

I_1_0 = pi^2/12;

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^3 / 3 + x(i)*pi^2/3 + I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % Для страховки добавим 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)));
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
   
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)));
end

grid on, hold on
plot(x, delta, 'k', 'linewidth', 2)
plot(x_minus, delta_minus, 'k', 'linewidth', 2)

k = 3; 

a = [1 ...
     0.1583482145380455955096383 ...
     0.0460645149909308107878344 ...
     0.0048861379108841469134267 ...
     0.0004336733305971515517559 ...
     0.0000173435613795895152436];
 
b = [1 ...
     0.0125148812047107612191739 ...
     0.0266693407000929631393759 ...
     0.0003285431094547362504004 ...
     0.0000820910787890062715299];

x_minus = 0:-0.1:-30;
y_minus = log(1 + exp(x_minus));

I_3_0 = 7*pi^4/120;

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^4 / 4 + x(i)*x(i)*pi^2/2 + 7*pi^4/60 - I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % Для страховки добавим 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)));
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
    
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)));
end

grid on, hold on
plot(x, delta, 'k', 'linewidth', 2)
plot(x_minus, delta_minus, 'k', 'linewidth', 2)

%% 

k = 4; 

a = [1 ...
     0.0560148791230902149024568 ...
     0.0351117957891800867706741 ...
     0.0021834386943672331415760 ...
     0.0002464861525522946634693 ...
     0.0000092228177886669241259];
 
b = [1 ...
    -0.0611726208769112866900252 ...
     0.0279968542816146833953639 ...
    -0.0007512148294307540141223 ...
     0.0000860680747142919882956];


x_minus = 0:-0.1:-30;
y_minus = log(1 + exp(x_minus));

% I_approximate = [];
a1 = (gamma(k+2))^(-1/k);
for i = 1:length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    I_precesion(i) = x(i)^5 / 5 + 2*x(i)*x(i)*x(i)*pi^2/3 + 7*x(i)*pi^4/15 + I_precesion_minus(i);
    
    I_approximate(i) = gamma(k+1)*y(i).*(1+a1*y(i))^k;
    % Для страховки добавим 10^-17
    delta(i) = log10(abs(I_approximate(i)/I_precesion(i)));
end

x_minus = -10:0.1:0;
y_minus = log(1 + exp(x_minus));

I_precesion_minus = [];

for i = 1:length(y_minus)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m))).^k;
    
    I_approximate_minus(i) = gamma(k+1)*y_minus(i).*(1+a1*y_minus(i))^k;
    % Для страховки добавим 10^-17
    delta_minus(i) = log10(abs(I_approximate_minus(i)/I_precesion_minus(i)));
end

grid on, hold on
xlabel('x'); ylabel('lg(d)');
plot(x, delta, 'k', 'linewidth', 2)
plot(x_minus, delta_minus, 'k', 'linewidth', 2)
text(11.05,0.07,' k = 1');
text(11.05,0.16,' k = 2');
text(11.05,0.26,' k = 3');
text(11.05,0.37,' k = 4');