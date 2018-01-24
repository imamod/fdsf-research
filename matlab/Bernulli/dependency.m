clc
clear all
close all

format long

n = 1:6;
dzeta = [pi^2/6 pi^4/90 pi^6/945 pi^8/9450 pi^10 / 93555 691*pi^12 / 638512875 ];
disp('Dzeta:')
disp(dzeta')
%%  Асимптотика чисел Бернулли и дзета=функции Римана
% B = [1/6 -1/30 1/42 -1/30 5/66 -691/2730 7/6 -3617/510 43867/798 -174611/330 ... 
%     854513/138 -236264091/2730 8553103/6 -23749461029/870 8615841276006/14322 ...
%     -7709321041217/510 2577867858367/6];
B = [1/6 -1/30 1/42 -1/30 5/66 -691/2730 7/6 -3617/510 43867/798 -174611/330 ... 
    854513/138]; %-236264091/2730 8553103/6 -23749461029/870 8615841276006/14322 ...
    %-7709321041217/510 2577867858367/6];
m = 1:length(B);
disp('from tables')
% disp(B')
% dzeta_test_stirling = (pi*exp(1)./m).^(2*m).*abs(B)./(4*sqrt(pi*m));
% disp('Stirling with 1 term:')
% disp(dzeta_test_stirling')

dzeta_check = ((-1).^(m-1)).*(2.^(2*m-1)).*(pi.^(2*m)).*B./gamma(2*m+1);
% disp('Dzeta via Bernulli:')
% disp(dzeta_check')

multiplier = ((-1).^(m-1)).*gamma(2*m+1)./((pi.^(2*m)).*(2.^(2*m-1)));
% B_test_as_article = multiplier./(1/2 + sqrt(1/4-((6/(pi^2))*(1-6/(pi^2))).^m));
B_test_as_article = multiplier.*(1 + 1./(2.^(2.*m)));
disp((B./B_test_as_article)')
disp('------------------------')

B_test = multiplier./(1 - 1./(2.^(2.*m)));

disp((B./B_test)')
disp('------------------------')
% result = (B./multiplier - 1 - 1./(2.^2.*m));
B_test = multiplier.*(1 + 2.^(-2*m)+3.^(-2*m)+4.^(-2*m)+5.^(-2*m));
result = B./B_test;
disp(result')
B_test = multiplier.*(1 + 2.^(-2*m)+3.^(-2*m)+4.^(-2*m)+5.^(-2*m) + 1./((5.5.^(2.*m-1)).*(2*m-1)));
result = B./B_test;
disp(result')

% disp('------------------------')
% result = (multiplier./B - 1 + 1./(2.^2.*m));
% disp(result')
% disp('------------------------')
% i = 1:length(m)-1;
% disp(((result(i)-1)./(result(i+1)-1))')
% disp('from assimptotic')
% disp(B_test')
% dzeta_test = ((-1).^(m-1)).*(2.^(2*m-1)).*(pi.^(2*m)).*B_test./gamma(2*m+1);
% disp('Dzeta via asimptotic Bernulli:')
% disp(dzeta_test')
% grid on, hold on
% plot(m, log10(dzeta_test-1), 'k*')
% disp((log10(dzeta_test-1)+2*m.*log10(2))')
% plot(m, log10(dzeta_test-1)+2*m.*log10(2), 'r*')
% disp((log10(dzeta_test-1)+m.*log10(2))')
% disp((dzeta_test./ (1+1./(2.^(2*m))))')


