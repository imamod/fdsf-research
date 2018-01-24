clc 
clear all 
close all
format long

k = [-1.5 -0.5, 0.5, 1.5, 2.5, 3.5];
x = 3:0.1:29.9;
y = log(1+exp(x));
C1 = k.*(k+1)*pi^2/6;

grid on, hold on
% while ( true )
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k-1.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (-2*I_transp*(-0.5)./(y.^(-0.5)));
%     z = (-2*I_transp.*(-0.5)./y.^(-0.5)).^(1/(-0.5));
    z = ((-2*I_transp.*(k(1)+1)./y.^(k(1) + 1)).^(2/k(1)) - 1).*(y.^2)*k(1)/(2*C1(1));
    plot(1./(y.^2), z, 'y*');
    disp('k = -3/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    disp('d2I');
    delta = ones(1, size(dI)-1);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_km32.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k-0.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = I_transp*(k(1)+1)./(y.^(k(1)+1));
%     z = (I_transp.*(k(1)+1)./y.^(k(1) + 1)).^(1/k(1));
    z = ((I_transp.*(k(2)+1)./y.^(k(2) + 1)).^(2/k(2)) - 1).*(y.^2)*k(2)/(2*C1(2));
    plot(1./(y.^2), z, 'b*')
%     count = count + 1;
% end
    disp('k = -1/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_km12.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k0.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp.*(k(2)+1)./(y.^(k(2)+1)));
%     z = (I_transp.*(k(2)+1)./y.^(k(2)+1)).^(1/k(2));
    z = ((I_transp.*(k(3)+1)./y.^(k(3) + 1)).^(2/k(3)) - 1).*(y.^2)*k(3)/(2*C1(3));
    plot(1./(y.^2), z, 'g*')
    disp('k = 1/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_k12.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k1.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp.*(k(3)+1)./(y.^(k(3)+1)));
%     z = (I_transp.*(k(3)+1)./y.^(k(3)+1)).^(1/k(3));
    z = ((I_transp.*(k(4)+1)./y.^(k(4) + 1)).^(2/k(4)) - 1).*(y.^2)*k(4)/(2*C1(4));
    plot(1./(y.^2), z, 'r*')
   
    disp('k = 3/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_k32.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k2.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp.*(k(4)+1)./(y.^(k(4)+1)));
%     z = (I_transp.*(k(4)+1)./y.^(k(4)+1)).^(1/k(4));
    z = ((I_transp.*(k(5)+1)./y.^(k(5) + 1)).^(2/k(5)) - 1).*(y.^2)*k(5)/(2*C1(5));
    plot(1./(y.^2), z, 'm*')
    disp('k = 5/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_k52.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_right_k3.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp.*(k(5)+1)./(y.^(k(5)+1)));
%     z = (I_transp.*(k(5)+1)./y.^(k(5)+1)).^(1/k(5));
    z = ((I_transp.*(k(6)+1)./y.^(k(6) + 1)).^(2/k(6)) - 1).*(y.^2)*k(6)/(2*C1(6));
    plot(1./(y.^2), z, 'k*')
    disp('k = 7/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(1/(y(i)^2)-1/(y(i-1)^2))];
    end
    disp('dI');
    disp(dI);
    for i=2:length(dI)
        delta(i-1) = dI(i)-dI(i-1);
    end
    f = fopen('d2I_k72.txt','w');
    fprintf(f, '%2.15f\n', delta);
    fclose(f);
    
    legend('-3/2', '-1/2', '1/2', '3/2', '5/2', '7/2');