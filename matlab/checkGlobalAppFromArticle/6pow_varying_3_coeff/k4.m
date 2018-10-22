clc
clear all
close all

format long

k = 4; N = 4;
n = 0:N+1; m = 0:N;

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
 
x = 0:0.1:50;
y = log(1 + exp(x));

x_minus = 0:-0.1:-50;
y_minus = log(1 + exp(x_minus));

c3 = ((k+1)*pi^2)/(gamma(k+2)^(6/k));
c4 = 1/(gamma(k+2)^(6/k));
%% c1 = 0.63, c2 = 0.44, c3 = 0.04
x_probe = [0.05; 2.3; 5.8];
y_probe = log(1+exp(x_probe));
if (x_probe(1) ~=0 || x_probe(2) ~= 0)
    y_probe_m = log(1+exp(-x_probe));
else
    y_probe_m = y_probe;
end
y_interp = [y_probe y_probe.^2 y_probe.^4];

if x_probe(1) ~=0 || x_probe(2) ~= 0
    P_probe = x_probe.^5/5 + 2*(x_probe.^3)*(pi^2)/3 + 7*x_probe.*(pi^4)/15;
end

for i = 1:length(y_probe)
    I_minus(i) = gamma(k+1)*y_probe_m(i).*(sum(a.*(y_probe_m(i).^n))/sum(b.*(y_probe_m(i).^m)))^k;
end
if x_probe(1) ~=0 || x_probe(2) ~= 0
    I = P_probe + ((-1)^k)*I_minus';
else 
end

Right = (I./(gamma(k+1)*y_probe)).^(6/k) - 1 - c4*y_probe.^6;
C = inv(y_interp)*Right;
c1 = C(1); c2 = C(2); c3 = C(3); %% ñ1 = 0.633462171838703; c2 = 0.440636818087044; c3 = 0.039244237631519
% c1 = 0.6334, c2 = 0.4406, c3 = 0.0392
P = x.^5/5 + 2*(x.^3)*pi^2/3 + 7*x.*pi^4/15;
for i = 1: length(y)
    I_precesion_minus(i) = gamma(k+1)*y_minus(i).*(sum(a.*(y_minus(i).^n))/sum(b.*(y_minus(i).^m)))^k;
    I_precesion(i) = P(i) + I_precesion_minus(i);

end

I_approximate = (1 + c1*y+c2*y.^2+c3*y.^4 + c4*y.^6);
I_approximate_minus = (1 + c1*y_minus+c2*y_minus.^2+c3*y_minus.^4 + c4*y_minus.^6);

% delta = (I_approximate./I_precesion) - 1;
% delta_minus = (I_approximate_minus./I_precesion_minus) - 1; 

delta = (I_approximate.*(gamma(k+1)*y./I_precesion).^(6/k) - 1)*(k/6);
delta_minus = (I_approximate_minus.*(gamma(k+1)*y_minus./I_precesion_minus).^(6/k) - 1)*(k/6);    
d = [delta_minus, delta];
hold on, grid on
plot([x_minus, x], d, 'k', 'linewidth',2)

% plot([-10,0;30,0],[0,-0.01;0,0.01],'k','linewidth',1.5);
% axis([-30 30 -2.5 1]),grid on
