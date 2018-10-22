clear all
close all
clc

k = 1:5;
dzeta = [pi^2/6 pi^4/90 pi^6/945 pi^8/9450 pi^10/93555 691*pi^12/638512875];

res = [];
for i = 2:length(dzeta)
    res = [res, (dzeta(i-1)-1)/(dzeta(i)-1)];
end
asimpFromDzeta = 2*((-1).^(k+1)).*(res.*gamma(2*k+1)./((2*pi).^(2*k)));
disp('Asimpfrom dzeta:');
disp(asimpFromDzeta)
k = 1:6;
% bernulli = [1 -1/2 1/6 0 -1/30 0 1/42 0 -1/30 0 5/66 0 -691/2730 0 7/6 0 -3617/510 0]'
disp('bernulli from wiki:')
bernulli = [1  1/6 -1/30 1/42 -1/30 5/66 -691/2730 7/6 -3617/510]';
disp(bernulli)

bernulli = 2*((-1).^(k+1)).*(dzeta.*gamma(2*k+1)./((2*pi).^(2*k)));
disp('bernulli:')
disp(bernulli')

res = [];
for i = 2:length(bernulli)
    res = [res, (bernulli(i-1))/(bernulli(i))];
end

disp(res')