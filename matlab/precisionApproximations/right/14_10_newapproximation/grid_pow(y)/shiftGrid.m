%%
%  * Замена переменных. Исследование сетки с автоматическим выравниванием экстремумов.
%  * ksi = 1/y; delta - экстремумы, полученные на сетке с ksi.
%  * eta(i) = ksi(i) + 0.5*(ksi(i+1)-ksi(i-1))*tau*(1+sqrt(delta(i-0.5)/delta(i+0.5)))/(1-sqrt(delta(i-0.5)/delta(i+0.5)))
%  * tau = 0.5

%   tau = 0.5;
%   tau= 1;
%   tau = 0.75;
function eta = shiftGrid(y, delta, tau)
    ksi = 1./y; 
    eta = [ksi(1)];
    for i = 2:length(ksi) - 1
        distance = (ksi(i + 1) - ksi(i - 1)) / 2;
        num = 1 - sqrt((delta(i - 1) / delta(i)));
        denom = (1 + sqrt(delta(i - 1) / delta(i)));
        eta = [eta, ksi(i) + tau*distance*num / denom];
    end
    eta = [eta, ksi(end)];
end
