function [TRUSS]=Truss_modal_analysis_37bar(Individual)
%PFJA Individual=[0.9724,	1.3603,	1.4976,	1.614,	1.6929,	0.00032214,	0.0001,	0.0001,	0.00026189,	0.0001389,	0.00015061,	0.00027945,	0.00015383,	0.00017091,	0.000264,	0.00011418,	0.00013366,	0.0002532,	0.0001];
%MJA Individual=[0.8912,	1.2439,	1.4101,	1.5595,	1.6312,	0.00036349,0.0001139,	0.0001,	0.00027251,	0.00013238,	0.0001349,	0.00030344,	0.00013394,	0.00017539,	0.00026325,	0.00012829,	0.00013327,	0.00030346,	0.0001];
%HSPO Individual=[0.9606,1.3425,1.5219,1.6567,1.733,0.00030179,0.0001,0.00010001,0.00025478,0.00012429,0.00012679,0.00025675,0.00014142,0.00015449,0.00025457,0.00012148,0.00013371,0.00023914,0.0001];
%AHEFA Individual=[0.9589,1.345,1.5355,1.6668,1.7397,0.0002821,0.00010019,0.00010001,0.00025308,0.0001221,0.00012429,0.00024718,0.00014018,0.00015061,0.00025604,0.00012146,0.00013605,0.00023992,0.0001];
%EACCS Individual=[0.9352,1.3114,1.4983,1.6284,1.7064,0.00029045,0.00010013,0.00010156,0.00026157,0.00011372,0.00012058,0.00025863,0.0001373,0.00014763,0.00024519,0.00012508,0.00013155,0.00022816,0.00010003];
%MS-DE Individual=[0.9522,1.3435,1.5387,1.675,1.7501,0.00028279,0.0001028,0.00010104,0.00023892,0.00012008,0.00012553,0.00024195,0.00013621,0.00015691,0.0002497,0.00012206,0.00014111,0.0002522,0.00010013];
%MJA  Individual=[0.8912,1.2439,1.4101,1.5595,1.6312,0.00036349,0.0001139,0.0001,0.00027251,0.00013238,0.0001349,0.00030344,0.00013394,0.00017539,0.00026325,0.00012829,0.00013327,0.00030346,0.0001];
%PFJA Individual=[0.9724,1.3603,1.4976,1.614,1.6929,0.00032214,0.0001,0.0001,0.00026189,0.0001389,0.00015061,0.00027945,0.00015383,0.00017091,0.000264,0.00011418,0.00013366,0.0002532,0.0001];

%AO [20,40,60] Individual=[4.02039016825054,2.26511868907031,3.70000417998373,3.98223389410482,2.50000152810858,0.0001,0.000100000342474174,0.000167206129746094,0.000180346967353315,0.000134047865367433,0.0001,0.000157488791003177,0.000165713781198278];

%SMA [20,40,60] Individual=[5.2854,1.635,3.798,3.6108,2.5908,0.00010116,0.0002092,0.00017732,0.00011536,0.00014278,0.0001,0.000186,0.00015233];

%LSMA [20,40,60]  Individual=[5.727,2.1605,3.7099,3.9476,2.5,0.0001,0.00013956,0.00013564,0.00016596,0.00015107,0.0001,0.00015806,0.000142];

%ESMA [20,40,60] Individual=[5.53653270122305,2.50678317963573,3.70000008877358,4.0707085304262,2.5,0.000100000028357757,0.000100023181723625,0.00012391644422531,0.000155515358932352,0.000134479080580897,0.000100002054960058,0.000175453891153738,0.00014873560341596];

%AOSMA [20,40,60] Individual=[5.75468445555751,2.33903636636903,3.70006357136164,4.02368287986739,2.50000010060151,0.000100000009345347,0.000112424342309701,0.000138609567142769,0.000155838371068613,0.000155285129908316,0.000100000090235289,0.000210172971312067,0.000101730187944319];

%ASMA [20,40,60] Individual=[5.43042326879736,2.03194385801416,3.70025727501922,3.87883909045967,2.5,5.02125542084371E-05,0.000165344975514668,0.000137116026575359,0.000167234520638223,0.000179974880267473,0.000068413357561033,0.00017332744687279,0.000130010432396011];



Truss_Data=TrussData_37bar_modal;

Truss_Data.dnkoor([3 19],2)=Individual(1);
Truss_Data.dnkoor([5 17],2)=Individual(2);
Truss_Data.dnkoor([7 15],2)=Individual(3);
Truss_Data.dnkoor([9 13],2)=Individual(4);
Truss_Data.dnkoor(11,2)=Individual(5);

Truss_Data.A([1 27])=Individual(6);
Truss_Data.A([2 26])=Individual(7);
Truss_Data.A([3 24])=Individual(8);
Truss_Data.A([4 25])=Individual(9);
Truss_Data.A([5 23])=Individual(10);
Truss_Data.A([6 21])=Individual(11);
Truss_Data.A([7 22])=Individual(12);
Truss_Data.A([8 20])=Individual(13);
Truss_Data.A([9 18])=Individual(14);
Truss_Data.A([10 19])=Individual(15);
Truss_Data.A([11 17])=Individual(16);
Truss_Data.A([12 15])=Individual(17);
Truss_Data.A([13 16])=Individual(18);
Truss_Data.A(14)=Individual(19);

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

for i=1:3
g(i)=(Truss_Data.limfre(i)/Frekans(i))-1;
gg(i)=g(i)*(g(i)>0);
end

%% Object Function
ObjVal=0;
for i=1:size(Truss_Data.eldn,1)
ObjVal=ObjVal+elboy(i)*Truss_Data.A(i)*Truss_Data.GS;
end
%% Penalized Obj. Func.
PEN=10^8;
Z=ObjVal;
for k=1:length(gg)
     Z=Z+ PEN*gg(k);
end
TRUSS.PENALIZED=Z;