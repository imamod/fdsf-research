%%% ������� ������� �� ������� U=0
clc
clear all
close all

% ����� �����
N = 15;

iteration_number_A_01= 22;
iterations_A_01 = 1:iteration_number_A_01;

iteration_number = 19;
iterations_A_optimal = 1:iteration_number;
%% ����������� tau �� ������ �������� ��� N = 15
tau_list_A_01 = [0.206 0.210 0.218 0.233 0.259 0.3 0.362 0.464 0.637 0.933];
tau_list_A_01 = [tau_list_A_01 ones(1, iteration_number_A_01 - length(tau_list_A_01))];

tau_list_optimal = [0.413, 0.438 0.537 0.815];
tau_list_optimal = [tau_list_optimal ones(1, iteration_number - length(tau_list_optimal))];
hold on
plot(iterations_A_01, tau_list_A_01, 'k.-', 'MarkerSize', 20);
plot(iterations_A_optimal, tau_list_optimal, 'k-o');
ylabel('\tau');xlabel('i');% ����� ��������
axis([1 iteration_number_A_01 0.2 1]);
%title('����������� \tau �� ������ ��������');
%% ����������� lg(L) �� ������ ��������

Lg_L_A_01 = [
   3.331993677923673
   3.144964877936638
   2.958505435134404
   2.765649070983328
   2.563260358762101
   2.347293908527719
   2.114123505544265
   1.853703636107110
   1.553494835883863
   1.194497390879759
   0.750761110052080
   0.373514371463697
   0.065362834747671
  -0.200092468018445
  -0.437817490309249
  -0.658600157736869
  -0.871075360469714
  -1.075546885385150
  -1.278472860911900
  -1.476210886507242
  -1.675776701760772
  -1.870323162296491];

Lg_L_optimal = [3.331993677923673
   2.952648925201254
   2.553252943301858
   2.098770991723886
   1.506456237615972
   0.959768446744644
   0.535343254148449
   0.198245505375607
  -0.087833108144852
  -0.334318941070932
  -0.562967389783021
  -0.777865740783973
  -0.985906220726220
  -1.189129479711414
  -1.389073691208328
  -1.587748088919476
  -1.784494537290194
  -1.981176560841204
  -2.176781181245253];
figure; hold on;
plot(iterations_A_01, Lg_L_A_01, 'k.-', 'MarkerSize', 20);
plot(iterations_A_optimal, Lg_L_optimal, 'k-o');
y_bottom = Lg_L_optimal(end)-0.1;
axis([1 iteration_number_A_01 y_bottom Lg_L_A_01(1)+0.1]);
xlabel('n'); ylabel('lg(L)');


%% ��������� �������� L(n)/L(n+1) - ����� �� ��������� (A = 0.1,0.2)
neighbours_A_01= [
   1.538256645519320
   1.536241322587496
   1.559036791442962
   1.593634460223865
   1.644244698448675
   1.710686400948768
   1.821460968986430
   1.996221828385346
   2.285585357019035
   2.778025838637857
   2.383673336785748
   2.033066277244877
   1.842702829738533
   1.728721456086606
   1.662580442444887
   1.631079775473323
   1.601295653862725
   1.595607156569342
   1.576659912585128
   1.583309491514114
   1.565115745313845
   1.580208101065270];

neighbours_A_optimal = [
   2.395216378531628
   2.508395320424986
   2.847619456266976
   3.911242596669263
   3.521176468524227
   2.657205805335140
   2.173190254330285
   1.932318062602226
   1.763948221971282
   1.692966828342691
   1.640205829123076
   1.614509035797583
   1.596699756224576
   1.584689613683019
   1.580062979249958
   1.573064203370970
   1.572830866625075
   1.568933804207316];
%figure; hold on;
%plot(iterations_A_01, neighbours_A_01, 'k.-', 'MarkerSize', 20);
%plot(iterations_A_optimal(2:end), neighbours_A_optimal, 'k-o');
%xlabel('n'); ylabel('L(n)/L(n+1)');
%% ������� �����������, ��������� � �������� beta*asinh(alpha*u)

clear first_it_profile; load('first_it_profile.mat', 'first_it_profile');
clear last_it_profile; load('last_it_profile.mat', 'last_it_profile');
clear first_it_ekstr; load('first_it_ekstr.mat', 'first_it_ekstr');
clear first_it_X; load('first_it_X.mat', 'first_it_X');
clear last_it_X; load('last_it_X.mat', 'last_it_X');
disp(first_it_ekstr);

alpha = 1/min(abs(first_it_ekstr));
beta = 1;

delta_first = beta.*asinh(alpha.*first_it_profile);
delta_final = beta.*asinh(alpha.*last_it_profile);

figure; hold on;
xlabel('x'); ylabel('\sigma(u)');
plot(first_it_X, delta_first, 'k-', 'MarkerSize', 20);
plot(last_it_X, delta_final, 'k-','linewidth', 2);
axis([0 1 -8.5 6.1]);
line([0,0;1,0],[0,-8.5;0,6.1], 'color','black', 'linewidth', 1.2);
%% �������������� ���������������
% ��� ������� N lg(min|u_ekstr|),lg(max|u_ekstr|) � lg(u_ekstr_final)

min_ekstr_it = [
    0.111111111111111   
    0.044326669921875    
    0.0143999999999997   
    0.00564106025270035  
    0.00191246844427052 
    0.00073928385972974  
    0.000256117068707873 
    9.82327499998314e-05 
    3.44299796043554e-05 
    1.31291358177561e-05 
    4.6379788053073e-06  
    1.76078534861826e-06 
    6.25544291531682e-07]';

max_ekstr_it = [
    0.197214814814815  
    0.113462958984375   
    0.0691539488640001  
    0.0436875444444448 
    0.0284273820236712
    0.0188032111889762 
    0.0125940687770854
    0.00851908369607181
    0.0058089071730371
    0.00398899008243158
    0.00276409142284511 
    0.0019238512695946
    0.00134416907829661]';

final_ekstr_it = [
    0.17211415895692
    0.0803314306019185
    0.0386288655241675
    0.0186936863044457
    0.00913405652458803
    0.00449357693453892
    0.00221317632959033
    0.00109621087053446
    0.000542254274531168
    0.000269179898440033
    0.00013355584571546
    6.6453130694289e-05
    3.3031008103579e-05]';

clear n; n = 4:N+1;
figure; hold on;
plot(n, log10(min_ekstr_it), 'k.-', 'MarkerSize', 20);
plot(n, log10(max_ekstr_it), 'k.-', 'MarkerSize', 20);
plot(n, log10(final_ekstr_it), 'ko-');
xlabel('N+1'); ylabel('lg(|p_n_-_1_/_2|)');

%% ����������� �� ��������� a = A � b = tau
%% ������ 2� ������ ���� �� � ����� DONE
clear n; n = 0:N;
grid_init = (2.*n-N)/N;
grid_final = [-1.000000000000000
  -0.961510406769359
  -0.886106889635691
  -0.776660170393637
  -0.637224576932842
  -0.473441980548909
  -0.291518946248502
  -0.098429017261113
   0.098429017261917
   0.291518946249287
   0.473441980548344
   0.637224576932356
   0.776660170392936
   0.886106889634653
   0.961510406770045
   1.000000000000000];
grid_size = length(grid_init);
y_value = zeros(1, grid_size);

figure; hold on;
%subplot(2,1,1);
plot(grid_init, y_value, 'k.-', 'MarkerSize', 20)
%axis([-1 1 -0.5 +0.5]);

%subplot(2,1,2);
plot(grid_final, y_value-0.001, 'ko-', 'MarkerSize', 8)
axis([0 1 -0.0015 +0.0005]);
