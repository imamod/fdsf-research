clear all
close all
clc

coefs = [pi^2/6 pi^4/90 pi^6/945 pi^8/9450 pi^10/93555 691*pi^12/638512875]

res = [];
for i = 2:length(coefs)
    res = [res, (coefs(i-1)-1)/(coefs(i)-1)];
end

disp(res')