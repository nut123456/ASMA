function Truss_Data=TrussData_37bar_modal

% Modulus of Elasticity
E=2.10e11;

% Node coordinates
dnkoor=zeros(20,3);
dnkoor=[0	0	0;
        1	0	0;
        1	0	0;
        2	0	0;
        2	0	0;
        3	0	0;
        3	0	0;
        4	0	0;
        4	0	0;
        5	0	0;
        5	0	0;
        6	0	0;
        6	0	0;
        7	0	0;
        7	0	0;
        8	0	0;
        8	0	0;
        9	0	0;
        9	0	0;
        10	0	0];

% Bar connections between nodes 
eldn=[
1	3	;
2	3	;
3	4	;
3	5	;
4	5	;
5	6	;
5	7	;
6	7	;
7	8	;
7	9	;
8	9	;
9	10	;
9	11	;
10	11	;
10	13	;
11	13	;
12	13	;
12	15	;
13	15	;
14	15	;
14	17	;
15	17	;
16	17	;
16	19	;
17	19	;
18	19	;
19	20	;
1	2	;
2	4	;
4	6	;
6	8	;
8	10	;
10	12	;
12	14	;
14	16	;
16	18	;
18	20	];

% Support conditions
MesKos=[1	1	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        0	0	1	;
        1	1	1	];

% Fixed cross sections
A=zeros(size(eldn,1),1);
A(28:37)=4e-3; %m2

% Frequency limits (Hz)
limfre(1)=20;
limfre(2)=40;
limfre(3)=60;

% Unit weight
GS=7800;%kg/m^3

% added masses
add_mass=zeros(size(dnkoor,1));
add_mass(2:2:18)=10; %kg

Truss_Data=struct('dnkoor',dnkoor,'MesKos',MesKos,'eldn',eldn,'E',E,'GS',GS,'add_mass',add_mass,'limfre',limfre,'A',A);