function max_delta = readData(k, N)
baseSize = 2*N+1;
% x_star = 3;
x_star = 5;
y_star = log(1+exp(x_star));

% задаем лин-триг сетку 
j = 1:baseSize;
alpha = 2/(2+pi);
y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));

prefix = strcat('linTrig/k12/n', num2str(N));
filename = strcat(prefix, '/Y.txt');
f = fopen(filename,'r');
Y = fscanf(f, '%f');
Y = Y';
fclose(f); 

filename = strcat(prefix, '/I_base.txt');
f = fopen(filename,'r');
I_base = fscanf(f, '%f');
I_base = I_base';
fclose(f);

filename = strcat(prefix, '/I_add.txt');
f = fopen(filename,'r');
I_add = fscanf(f, '%f');
I_add = I_add';
fclose(f);

max_delta = commonPart(y0, Y, I_base, I_add, k, N);
end
