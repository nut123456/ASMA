function [TRUSS]=Truss_modal_analysis_52bar(Individual)
%Individual=[5.9463,2.2264,3.7165,3.9496,2.5,1.0063/10000,1.2140/10000,1.2734/10000,1.4194/10000,1.4147/10000,1/10000,1.5435/10000,1.3919/10000];%GWO
%SMA% Individual=[5.7036,	2.1294,	3.7268,	3.9222,	2.5001,	0.0001,
%0.00013971,	0.00012738,	0.00015981,	0.00014042,	0.0001,	0.00016354,	0.00015141]

%LSMA% Individual=[5.8221,	2.1102,3.7604	,3.9616	,2.5,	0.00010211,0.0001573,	0.00014216,	0.00013456,	0.00014141,	0.00010029,	0.00015474,	0.00015685]

%CELSMA Individual=[4.0063,	2.213,	3.7001,	3.9448,	2.5,	0.00010016,0.00010009,	0.00014864,	0.00019687,	0.00012841,	0.00010011,	0.00015517,	0.00018154]
%Individual=	[5.7514,	2.2353,	3.7198,	3.9534,	2.4944,	0.00005623,0.0001167,	0.00011772,	0.00014399,	0.00015879,	0.000082,	0.0001628,	0.0001363];
%ASMA199.2923 Individual=[5.809014087	,2.293320324	,3.700330651	,4.015362201,2.500067349	,7.64E-05	,0.000129607	,1.35E-04	,0.000151977	,0.000182277	,6.97E-05	,0.000181644	,1.20E-04	];
%ASMA198.4584134Individual=[5.38452064986677,2.00211113292283,3.70188780537047,3.84058974085702,2.50001072095218,4.34543469004913E-05,0.000189310886479186,0.000128318650315912,0.000157656886848597,0.000168981180851023,6.64580397743748E-05,0.000167411060318366,0.000153247097822239]
%ASMA197.8171Individual=[5.34495940816609,2.08446460274263,3.72128278247421,3.85795992526574,2.52654328138654,5.34735670927962E-05,0.000158100228009135,0.000132860956419954,0.000159105651119324,0.000158614724956968,7.12371298337564E-05,0.000171975721404467,0.00014840622152288];
%AO [15.9078,28.648228.6482,28.7235,29.7542] Individual=[4.02039016825054,2.26511868907031,3.70000417998373,3.98223389410482,2.50000152810858,0.0001,0.000100000342474174,0.000167206129746094,0.000180346967353315,0.000134047865367433,0.0001,0.000157488791003177,0.000165713781198278];

%SMA [    9.5794   28.6671   28.6852   28.6852   29.3948] Individual=[5.2854,1.635,3.798,3.6108,2.5908,0.00010116,0.0002092,0.00017732,0.00011536,0.00014278,0.0001,0.000186,0.00015233];

%LSMA [10.8540   28.6478   28.6478   28.6482   28.9313]  Individual=[5.727,2.1605,3.7099,3.9476,2.5,0.0001,0.00013956,0.00013564,0.00016596,0.00015107,0.0001,0.00015806,0.000142];

%ESMA [ 15.0763   28.6486   28.6489   28.6489   29.6985] Individual=[5.53653270122305,2.50678317963573,3.70000008877358,4.0707085304262,2.5,0.000100000028357757,0.000100023181723625,0.00012391644422531,0.000155515358932352,0.000134479080580897,0.000100002054960058,0.000175453891153738,0.00014873560341596];

%AOSMA [   11.7957   28.6481   28.6481   28.6494   29.1056] Individual=[5.75468445555751,2.33903636636903,3.70006357136164,4.02368287986739,2.50000010060151,0.000100000009345347,0.000112424342309701,0.000138609567142769,0.000155838371068613,0.000155285129908316,0.000100000090235289,0.000210172971312067,0.000101730187944319];

%ASMA [ 10.1140   28.6480   28.6480   28.6550   28.8525] Individual=[5.43042326879736,2.03194385801416,3.70025727501922,3.87883909045967,2.5,5.02125542084371E-05,0.000165344975514668,0.000137116026575359,0.000167234520638223,0.000179974880267473,0.000068413357561033,0.00017332744687279,0.000130010432396011];

%MPA  [ 10.1328   27.8588   28.0840   28.0840   28.6421 ]199.4440 Individual=[5.2232,1.7915,3.7001,3.6571,2.5263,0.00010144,0.00015378,0.00012636,0.00016659,0.00012698,0.00010065,0.00016682,0.0001593,0.019944];
%IMPA Individual=[5.2872,2.3564,3.7119,3.9948,2.5,0.0001,0.00010002,0.00012728,0.00015005,0.00013597,0.0001,0.0001591,0.00014446,0.019483];




Truss_Data=TrussData_52bar_modal;

% ZA=Individual(5)+Individual(3)+Individual(1);
% XB=Individual(2);
% ZB=Individual(5)+Individual(3);
% XF=Individual(2)+Individual(4);
% ZF=Individual(5);

ZA=Individual(1);
XB=Individual(2);
ZB=Individual(3);
XF=Individual(4);
ZF=Individual(5);

Truss_Data.dnkoor=[0	0	ZA;
XB	0	ZB;
0	XB	ZB;
-XB	0	ZB;
0	-XB	ZB;
XF	0	ZF;
XF*cos(45/180*pi)	XF*cos(45/180*pi)	ZF;
0	XF	ZF;
-XF*cos(45/180*pi)	XF*cos(45/180*pi)	ZF;
-XF	0	ZF;
-XF*cos(45/180*pi)	-XF*cos(45/180*pi)	ZF;
0	-XF	ZF;
XF*cos(45/180*pi)	-XF*cos(45/180*pi)	ZF;
6	0	0;
6*cos(45/180*pi)	6*cos(45/180*pi)	0;
0	6	0;
-6*cos(45/180*pi)	6*cos(45/180*pi)	0;
-6	0	0;
-6*cos(45/180*pi)	-6*cos(45/180*pi)	0;
0	-6	0;
6*cos(45/180*pi) -6*cos(45/180*pi) 0];

Truss_Data.A(1:4)=Individual(6);
Truss_Data.A(5:8)=Individual(7);
Truss_Data.A(9:16)=Individual(8);
Truss_Data.A(17:20)=Individual(9);
Truss_Data.A(21:28)=Individual(10);
Truss_Data.A(29:36)=Individual(11);
Truss_Data.A(37:44)=Individual(12);
Truss_Data.A(45:52)=Individual(13);

if ZB==ZF && XB==XF
    ObjVal=10^8;
    gg(1:2)=10^8;
    
else 
    
dnsay=size(Truss_Data.dnkoor,1);
elsay=size(Truss_Data.eldn,1);

topser=dnsay*3;
yer=zeros(topser,1);

[rijitlik,elboy]=stiffness_3D_truss(topser,elsay,Truss_Data.eldn,Truss_Data.dnkoor,Truss_Data.E,Truss_Data.A);
[mass_mat]=mass_truss(topser,elsay,Truss_Data.eldn,Truss_Data.dnkoor,Truss_Data.GS,Truss_Data.A);

for i=1:dnsay
mass_mat(3*i-2,3*i-2)=mass_mat(3*i-2,3*i-2)+Truss_Data.add_mass(i);
mass_mat(3*i-1,3*i-1)=mass_mat(3*i-1,3*i-1)+Truss_Data.add_mass(i);
mass_mat(3*i,3*i)=mass_mat(3*i,3*i)+Truss_Data.add_mass(i);
end

tutser=find(Truss_Data.MesKos'==1);
ser=setdiff([1:topser]',[tutser]);
[FI,LAMDA]=eig(inv(mass_mat(ser,ser))*rijitlik(ser,ser));
[LAMDA_sorted, ind]=sort(diag(LAMDA),'ascend');


Frekans=sqrt(LAMDA_sorted)/(2*pi)

%% Constraints

g(1)=(Frekans(1)/Truss_Data.limfre(1))-1;
gg(1)=g(1)*(g(1)>0);

g(2)=(Truss_Data.limfre(2)/Frekans(2))-1;
gg(2)=g(2)*(g(2)>0);

%% Object Function
ObjVal=0;
for i=1:size(Truss_Data.eldn,1)
ObjVal=ObjVal+elboy(i)*Truss_Data.A(i)*Truss_Data.GS;
end

end
%% Penalized Obj. Func.
PEN=0;
Z=ObjVal;
for k=1:length(gg)
    out = imag(gg(k)) ~= 0;
    if out == true
        gg(k)=abs(gg(k));
    end
     Z=Z+ PEN*gg(k);
end
TRUSS.PENALIZED=Z;