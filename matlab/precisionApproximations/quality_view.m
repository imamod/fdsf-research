clc 
clear all 
close all
format long

% коэффициент k функции Ферми-Дирака
% k = -3/2;

% N = 6;
% x_star = 3;
% y_star = log(1+exp(x_star));

% a0 = 1; b0 = 1;
% a = zeros(1,N+1); b = zeros(1,N);

% задаем лин-триг сетку
% baseSize = 2*N+1;
% j = 1:baseSize;
% alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
% x0 = log(exp(y0)-1);

k = [-0.5, 0.5, 1.5, 2.5, 3.5];
x = -5:0.1:5;
y = log(1+exp(x));

% file_names = [
%     'D:\repositories\OpenFFD\testOpenFFD\I_k-0.5.txt'
%     'D:\repositories\OpenFFD\testOpenFFD\I_k0.5.txt'
%     'D:\repositories\OpenFFD\testOpenFFD\I_k1.5.txt'
%     'D:\repositories\OpenFFD\testOpenFFD\I_k2.5.txt'
%     'D:\repositories\OpenFFD\testOpenFFD\I_k3.5.txt'];
% I = size(:,file_names.size());
count = 1;
grid on, hold on
% while ( true )
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k-1.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (-2*I_transp./(gamma(-0.5)*y)).^(1/(-1.5));
%     z = (-2*I_transp./(gamma(-0.5)*y));
    z = log(-2*I_transp./(gamma(-0.5)*y));
    plot(y, z, 'y*');
    disp('k = -3/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k-0.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp./(gamma(k(1)+1)*y)).^(1/k(1));
%     z = I_transp./(gamma(k(1)+1)*y);
    z = log(I_transp./(gamma(k(1)+1)*y));
    plot(y, z, 'b*')
%     count = count + 1;
% end
    disp('k = -1/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k0.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp./(gamma(k(2)+1)*y)).^(1/k(2));
%     z = I_transp./(gamma(k(2)+1)*y);
    z = log(I_transp./(gamma(k(2)+1)*y));
    plot(y, z, 'g*')
    disp('k = 1/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k1.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp./(gamma(k(3)+1)*y)).^(1/k(3));
%     z = I_transp./(gamma(k(3)+1)*y);
    z = log(I_transp./(gamma(k(3)+1)*y));
    plot(y, z, 'r*')
   
    disp('k = 3/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end
    
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k2.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp./(gamma(k(4)+1)*y)).^(1/k(4));
%     z = I_transp./(gamma(k(4)+1)*y);
    z = log(I_transp./(gamma(k(4)+1)*y));
    plot(y, z, 'm*')
    disp('k = 5/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end
    f = fopen('D:\repositories\OpenFFD\testOpenFFD\I_k3.5.txt','r');
    I = fscanf(f, '%f');
    fclose(f);
    I_transp = transpose(I);
%     z = (I_transp./(gamma(k(5)+1)*y)).^(1/k(5));
%     z = I_transp./(gamma(k(5)+1)*y);
    z = log(I_transp./(gamma(k(5)+1)*y));
    plot(y, z, 'k*')
    disp('k = 7/2: ')
    dI=[];
    for i=2:length(y)
        dI = [dI, (z(i)-z(i-1))/(y(i)-y(i-1))];
    end
    for i=2:length(dI)
        disp(dI(i)-dI(i-1))
    end