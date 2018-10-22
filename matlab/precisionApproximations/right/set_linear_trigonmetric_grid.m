function [y0_inv, y0, x0] = set_linear_trigonmetric_grid(baseSize)
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);