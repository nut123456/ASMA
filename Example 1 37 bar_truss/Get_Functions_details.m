%% Benchmark Test functions
function [lb, ub, dim, fobj] = Get_Functions_details(F)
    switch F
        case 'F1'
            fobj = @F1;
            lb = [0.3, 0.3, 0.3, 0.3, 0.3, 0.5e-4, 0.5e-4,0.5e-4,0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4,  0.5e-4, 0.5e-4];
            ub = [3  , 3,   3,   3,   3  , 10e-4,  10e-4, 10e-4, 10e-4,   10e-4, 10e-4,  10e-4,  10e-4,  10e-4,  10e-4,   10e-4,  10e-4,   10e-4,   10e-4];

                %             LB=ones(1,19);
                % LB(1,1:5)=LB(1,1:5)*0.3;   %m
                % LB(1,6:19)=LB(1,6:19)*1e-4;  %m2
                % 
                % UB=ones(1,19);
                % UB(1,1:5)=UB(1,1:5)*3;   %m
                % UB(1,6:19)=UB(1,6:19)*10e-4;  %m2 

            dim = length(lb);
    end
end

function W = F1(Individual)

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


Frekans=sqrt(LAMDA_sorted)/(2*pi);

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
W = Z; % Assign value to W to ensure output is returned correctly
end