%% � ���� ����� �������� z � I ��������� ����������������, �.� 
% z = (Ik/(k!*y))^(1/k)
%
clc 
clear all 
close all
format long

% ����������� k ������� �����-������
k = 1/2;

N = 20;

x_star = 10;
y_star = log(1+exp(x_star));

y_star_inv = 1/y_star^2;

% ������ ���-���� �����
baseSize = N;
j = 1:baseSize;
alpha = 2/(2+pi);
% y0 = 0.5*y_star*(2*alpha*j/(baseSize)+(1-alpha)*(1 - cos(pi*j/(baseSize))));
y0_inv = 0.5*y_star_inv*(2*alpha*j/baseSize +(1-alpha)*(1 - cos(pi*j/baseSize)));
y0 = 1./sqrt(y0_inv(end:-1:1));
x0 = log(exp(y0)-1);

Y = [10.00004539889921062
10.01061241500150167
10.02121300045201302
10.03184733336678391
10.04251559318776899
10.05321796069554630
10.06395461802218527
10.07472574866427095
10.08553153749607389
10.09637217078287108
10.10724783619443734
10.11815872281869133
10.13261658948197486
10.14713661045560933
10.16171923235496166
10.17636490630122559
10.19107408798002901
10.20584723770099878
10.22068482045824922
10.23558730599186717
10.25055516885036155
10.26558888845413975
10.28068894916001241
10.29940692067547126
10.31822750490168872
10.33715164282356014
10.35618028755117059
10.37531440452141318
10.39455497170372489
10.41390297981003954
10.43335943250904485
10.45292534664486794
10.47260175246029412
10.49238969382460773
10.51583129822064500
10.53943072433993677
10.56318975107212843
10.58711018550511263
10.61119386350234528
10.63544265029467084
10.65985844108708491
10.68444316168090857
10.70919876911180779
10.73412725230414289
10.75923063274217206
10.78798241458730089
10.81696593452902455
10.84618432238670316
10.87564076748041764
10.90533852009327376
10.93528089297785755
10.96547126290844965
10.99591307228063997
11.02660983076005863
11.05756511698202971
11.08878258030400232
11.12359228941115497
11.15873189231667872
11.19420663286725315
11.23002187235191229
11.26618309290560305
11.30269590103407040
11.33956603126524421
11.37679934993248665
11.41440185909540794
11.45237970060417787
11.49073916031356823
11.53256819152834467
11.57485737685926175
11.61761521545363784
11.66085042786866310
11.70457196354289131
11.74878900857825670
11.79351099384791723
11.83874760344619581
11.88450878349781803
11.93080475134464535
11.97764600512911493
12.02774352892919474
12.07847496888048155
12.12985380767386978
12.18189393291286393
12.23460965288366609
12.28801571308252782
12.34212731354367421
12.39696012701407746
12.45253031802437960
12.50885456290873954
12.56595007082980686
12.62596082556740740
12.68683964843358858
12.74860767127404415
12.81128675321249055
12.87489951315201253
12.93946936406946158
13.00502054922094430
13.07157818038562347
13.13916827828479406
13.20781781532400601
13.27755476081769892
13.34968053320349313
13.42299461741884414
13.49753000616666299
13.57332098898817030
13.65040321854641370
13.72881378109810058
13.80859127147012444
13.88977587288512971
13.97240944201113777
14.05653559964413901
14.14219982746994475
14.22944957139289102
14.31833435196469040
14.40890588249726889
14.50121819549963575
14.59532777814108506
14.69129371751222202
14.78917785653262307
14.88904496143965517
14.99096290188917280
15.09500284480597365
15.20123946324206798
15.30783206277899033
15.41669888907995123
15.52792198006935820
15.64158757744821315
15.75778640776660922
15.87661398681215630
15.99817094963517050
16.12256340880209038
16.24990334377578094
16.38030902467014727
16.51390547402393594
16.64594922533396115
16.78121191930947020
16.91982650290994883
17.06193373957879444
17.20768281021167923
17.35723197158038644
17.51074927881883880
17.66841337947037971
17.83041438762677444
17.99695484788480471
18.16825080023482997
18.33495729749657244
18.50633859940518278
18.68261735744848195
18.86403135789912966
19.05083487076531057
19.24330014866226790
19.44171909570237489
19.64640512967617170
19.85769526455118239
20.07595244477190732
20.30156816815230059
20.51758378616541734
20.74064491306317493
20.97114306047913601
21.20950089137576811
21.45617548103599503
21.71166200522931788
21.97649792329980301
22.25126773663075852
22.53660841841941220
22.83321562963712026
23.14185085934676067
23.43218006246147667
23.73371789815274369
24.04720465328694345
24.37345092462940954
24.71334644434559280
25.06787029838216085
25.43810280521280021
25.82523938350977133
26.23060681468107092
26.65568240491510110
27.10211667818868264
27.51396064170514677
27.94516904978093663
28.39730856821410043
28.87212925585799184
29.37159312271353429
29.89790831847847130
30.45357031233302436
31.04141181889185219
31.66466375399796362
32.32703022089727796
33.03278151038990984
33.66997208372239214
34.34551840407664969
35.06342968795556914
35.82832750191778359
36.64557151028881066
37.52141829611779400
38.46322405306800363
39.47970622907895688
40.58128551912201232
41.78053909701657176
43.09281054260868871
44.25251010058249790
45.51117740881480955
46.88374025001036216
48.38848166362495817
50.04807630708657484
51.89104698825710926
53.95385876964260063
56.28400910358917031
58.94472710985684927
62.02237634004617206
65.63861499221097517
68.84236018546482683
72.56621922925350532
76.96809855311218485
82.28235842777056064
88.87510483730594046
97.35779944005794562
108.84932884388027219
125.68837861824874835
153.93619710622436969
217.69865768776054438
]';

X = log(exp(Y)-1); 

I_base = [21.34447149235519703
21.71749154180126595
22.23441576024407595
22.91397578719722006
23.78050192280538866
24.86581037374407899
26.21185815468908586
27.87458445885761549
29.92962375857960211
32.48105478190185380
35.67524537583092581
39.72361157264049325
44.94176410089801266
51.82065903192580691
61.16513411742646156
74.38867446557846108
94.22003233687448187
126.71195672139141664
188.71418053402453552
354.62734410341545299
]';

C1 = (k+1)*k*(pi^2)/6;
z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);

I_add = [21.34447149235516505
21.37775292147840744
21.41115804446531001
21.44468758395986541
21.47834226837277427
21.51212283193986252
21.54603001478113455
21.58006456296070752
21.61422722854735667
21.64851876967589916
21.68293995060922086
21.71749154180122687
21.76330458214054531
21.80934804630069834
21.85562375749185193
21.90213355864410616
21.94887931267860637
21.99586290278280032
22.04308623269067269
22.09055122696731033
22.13825983129865804
22.18621401278568683
22.23441576024406174
22.29421641721273417
22.35440063041730951
22.41497227321659125
22.47593527248294976
22.53729360954212524
22.59905132113326331
22.66121250038920110
22.72378129783815837
22.78676192242677345
22.85015864256566331
22.91397578719721295
22.98965533932097216
23.06593102348954361
23.14281023971070539
23.22030051377404192
23.29840949997099742
23.37714498388679019
23.45651488526554829
23.53652726095217673
23.61719030791216056
23.69851236633236269
23.78050192280538866
23.87452680448638986
23.96943834919976624
24.06524974616432999
24.16197445349139983
24.25962620516625279
24.35821901824962765
24.45776720030858797
24.55828535708442573
24.65978840040726183
24.76229155636598733
24.86581037374407188
24.98141568639174181
25.09830282567715187
25.21649423438646309
25.33601289438883697
25.45688234314294363
25.57912669081907708
25.70277063806414830
25.82783949443810201
25.95435919755199095
26.08235633293850242
26.21185815468908586
26.35332312853328318
26.49660882651263805
26.64175229313953963
26.78879160807005633
26.93776592301473727
27.08871550025462582
27.24168175284521709
27.39670728659551813
27.55383594391531332
27.71311284962850507
27.87458445885761549
28.04763514416358916
28.22324765671919167
28.40148202298053093
28.58240020332112508
28.76606617163155022
28.95254599891941893
29.14190794114834659
29.33422253156991033
29.52956267782099076
29.72800376407797884
29.92962375857961632
30.14203803995879838
30.35804701833908936
30.57774710776855898
30.80123828249804063
31.02862424514761486
31.26001260458481568
31.49551506418003299
31.73524762115656728
31.97933077781125988
32.22788976544377704
32.48105478190187512
32.74359913076091999
33.01120398579044490
33.28402416797493402
33.56222102865799428
33.84596280239398425
34.13542498314337337
34.43079072564592735
34.73225127397424217
35.04000641945189187
35.35426499032906378
35.67524537583095423
36.00317608744676789
36.33829636060542612
36.68085680019131445
37.03112007369948344
37.38936165620987140
37.75587063178839742
38.13095055639978170
38.51492038794857820
38.90811548966253497
39.31088871370341309
39.72361157264048614
40.13918019227890710
40.56512229502089895
41.00185035291170266
41.44979952752752439
41.90942927438545951
42.38122508675917999
42.86570039333956572
43.36339862591951544
43.87489547527078315
44.40080135564227959
44.94176410089800555
45.47860769771215672
46.03076409181691986
46.59893055594494626
47.18384838266968018
47.78630646464792164
48.40714523349043219
49.04726100017773405
49.70761074594389584
50.38921741951847366
51.09317580472833242
51.82065903192576428
52.53196893450618177
53.26662110143293205
54.02584157645008389
54.81094591927068649
55.62334764830645639
56.46456766682996431
57.33624480982585681
58.24014767124537428
59.17818789807605384
60.15243516948529390
61.16513411742639050
62.14005197078239462
63.15218902965055747
64.20382696504157138
65.29744252054834419
66.43572913001180780
67.62162150443272424
68.85832367861017644
70.14934110334442607
71.49851748587185796
72.91007722502408228
74.38867446557840424
75.78863164025672461
77.25187411992710906
78.78301519803837039
80.38713921410680996
82.06986418096957436
83.83741478180819229
85.69670781236142432
87.65545263523894448
89.72226984030059782
91.90683211068765956
94.22003233687446766
96.37101021039202919
98.64047778027121183
101.03902301985344536
103.57856735884483612
106.27258579646021985
109.13637258215523218
112.18736395174360609
115.44553286099065303
118.93387533313403992
122.67901443163154340
126.71195672139141664
130.39039386038382418
134.32849820507146887
138.55628146758328967
143.10869292793600493
148.02669563982161094
153.35864026184717090
159.16203802423106595
165.50587620953888290
172.47368198914136883
180.16763547586290883
188.71418053402453552
196.37637644495504219
204.80688679269391628
214.13413240716840846
224.51770092007890867
236.15859168181054883
249.31384085452771160
264.31790200748110919
281.61476593421929238
301.80775938347613874
325.73965757470205062
354.62734410341545299
380.89513603269688247
412.20403462433438335
450.26202874903549400
497.67746690538342591
558.65895414107694705
640.50333083812324730
757.16860197248990971
939.47524976382135264
1273.33442901271018854
2141.42896746656197138
]';

I = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

% ������ ������� A � B
B = z(1,:) - ones(1,baseSize);
A = zeros(baseSize,baseSize);
for i = 1:size(A,2)
    for j = 1:size(A,2)
        A(i,j) = y0(i)^(-2*j);
    end
end
disp('A matrix: ')
disp(A); 
disp('B matrix: ')
disp(B');
E = A\B';

a0 = 1;
a = zeros(1,baseSize);
for j = 1:length(E)
    a(j) = E(j);
end
disp('--------------------------------');
disp('������������ �:'); disp(a');
disp('----------------------------');

F_base = zeros(1,baseSize);
delta_base = zeros(1,baseSize);
for j = 1:baseSize
    F_base(j) = 1; 
    for n = 1:baseSize
        F_base(j) = F_base(j) + a(n).*y0(j).^(-2*n);
    end
    delta_base(j) = F_base(j)/z(j)-1;
end
%---------------------------------------
% ������� ��������������� �����

F = zeros(1,length(Y));
delta_additional = zeros(1,length(Y));
for j = 1:length(Y)
    F(j) = 1;
    for n = 1:baseSize
        F(j) = F(j) + a(n).*Y(j).^(-2*n);
    end
    delta_additional(j) = F(j)/I(j)-1;
end

disp('-------------------');
disp('����������:');
disp(max(abs(delta_additional(1:11))));
disp(max(abs(delta_additional(11:22))));
disp(max(abs(delta_additional(22:33))));
disp(max(abs(delta_additional(33:44))));
disp(max(abs(delta_additional(44:55))));
disp(max(abs(delta_additional(55:66))));
disp(max(abs(delta_additional(66:77))));
disp(max(abs(delta_additional(77:88))));
disp(max(abs(delta_additional(88:99))));
disp(max(abs(delta_additional(99:110))));
disp(max(abs(delta_additional(110:121))));
disp(max(abs(delta_additional(121:132))));
disp(max(abs(delta_additional(132:143))));
disp(max(abs(delta_additional(143:154))));
disp(max(abs(delta_additional(154:165))));
disp(max(abs(delta_additional(165:176))));
disp('-------------------');

grid on, hold on
xlabel('y'); %ylabel('d*10^1^0');
plot(Y,delta_additional, 'k','linewidth', 2.5);
plot(y0,delta_base, 'k*','linewidth',5)
% axis([0 y_star -1.5*10^(-1) 1.5*10^(-1)])
% line([0;y_star],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;log(1+exp(x_star))],[0; 0],'linewidth', 2, 'color', 'black');
% line([0;0],[-1.25; 1.35],'linewidth', 3, 'color', 'black');
% title({'�������-������������������ �����';['k = ', num2str(k), ', N = ', num2str(N), ', x_d_i_v = ', num2str(x_div)]})
% plot(log(1+exp(x_div)),2*10^-16,'b*');