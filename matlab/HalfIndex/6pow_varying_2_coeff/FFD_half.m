function I = FFD_half(ksi, x, k, a)
 
   exp_ksi = exp(-a.*ksi*ksi / (1 - ksi*ksi));
   
   I = (2*(a.^(k+1)).*(ksi^(2*k + 1)).*exp_ksi) ./ ...
            (((1 - ksi*ksi)^(k + 2)).*(exp_ksi + exp(-x)));
end