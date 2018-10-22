function max_delta = readData(k, N)
baseSize = 2*N;
% x_star = 2;
x_star = 3;
y_star = log(1+exp(x_star));

y_star_inv = 1/(y_star^2);

% задаем лин-триг сетку 
j = 1:baseSize;
alpha = 2/(2+pi);
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./sqrt(y0_inv(end:-1:1));
% x0 = log(exp(y0)-1);

prefix = strcat('k12/n', num2str(N));
filename = strcat(prefix, '/Y_12.txt');
f = fopen(filename,'r');
Y = fscanf(f, '%f');
Y = Y';
fclose(f); 

filename = strcat(prefix, '/I_base_12.txt');
f = fopen(filename,'r');
I_base = fscanf(f, '%f');
I_base = I_base';
fclose(f);

filename = strcat(prefix, '/I_add_12.txt');
f = fopen(filename,'r');
I_add = fscanf(f, '%f');
I_add = I_add';
fclose(f);

max_delta = commonPart(y0, Y, I_base, I_add, k, N);
end
