clc
clear all
close all
%%
N_trapz = [2 4 8 16 32 64 128 256];
I_trapz = [ 0.37130429265798032
            0.28472329634826737
            0.29051436415184395
            0.29050094109582197
            0.29050089612216462
            0.29050089616991714
            0.29050089616991764
            0.29050089616991757
          ];
      
I_trapz_new = [0.37130429265798032
0.28472329634826743
0.29051436415184406
0.29050094109582197
0.29050089612216456
0.29050089616991726
0.29050089616991759
0.29050089616991759];

I_prec = 0.29050089616991753;
figure, grid on, hold on
plot(N_trapz, log10(I_trapz./I_prec - 1),'k-o','MarkerSize', 7, 'linewidth', 2);
plot(N_trapz, log10(I_trapz_new./I_prec - 1),'r-o','MarkerSize', 7, 'linewidth', 2);
xlabel('N_m'); ylabel('ln |R_m|');
% title('trapz')
%%
N_quad = [2 4 8 16 32 64 128];
I_quad = [0.19814230003855449
          0.29630543195542058
          0.29048751803979989
          0.29050085114850727
          0.29050089621766995
          0.29050089616991781
          0.29050089616991759
         ];
     
% figure, grid on, hold on
% plot(N_quad,log10(I_quad./I_prec - 1),'k-s','MarkerSize', 7, 'linewidth', 2);
% xlabel('N_m'); ylabel('ln |R_m|');
% title('quad')
%%
N_simpson = [4 8 16 32 64 128 256];
I_simpson = [0.2558629642450298
            0.29244472008636951
            0.29049646674381457
            0.29050088113094552
            0.29050089618583486
            0.29050089616991759
            0.29050089616991759
             ];

% figure, grid on, hold on
% plot(N_simpson,log10(I_simpson./I_prec - 1),'k-d','MarkerSize', 7, 'linewidth', 2);
% xlabel('2N_m'); ylabel('ln |R_m|');
% title('simpson
% axis([0 64 -16 0]);