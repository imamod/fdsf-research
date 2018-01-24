function I = trapz(a, N, k, x)
     h = 1/N;
     u0 = FFD_half(0, x, k, a);
     I = u0 / 2; 

     for  i = 1:2:(N-1)/ 2
         I = I + FFD_half(i*h, x, k, a) + FFD_half((N - i)*h, x, k, a);
     end

     if (N == 2) 
         I = I + FFD_half(N*h/2, x, k, a);
     end

     I = h*I;
end