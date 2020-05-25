%%
%Module 1: Free stream conditions
function output=Moudgalya_Sandeep_Project(m2)

%Accepting all inputs
% z=input("Enter flight altitude: ");
% m1=input("Enter flight Mach no: ");
% nd=input("Enter inlet/diffuser efficiency: ");
% m2=input("Enter Mach no at diffuser exit: ");
% Tt3max=input("Enter maximum allowable total temperature at combustor exit: ");
% qf=input("Enter heating value of the fuel: ");
% nn=input("Enter nozzle efficiency: ");
% Ae=input("Enter nozzle exit area: ");
%Value for iteration
%z=4300;
%Values for validation:
% z=4300;
% m1=2.4;
nd=0.92;
% %m2=0.15;
% %m2=0.4;
Tt3max=2400;
qf=43.2*1000000;
nn=0.94;
Ae=0.015;
%Hypersonic Flight:
m1=5;
z=27400;

%Given values:
Ts=288;
ps=101.3*1000;
gamma=1.4;
zstar=8404;
R=286.9;
a=986;
b=0.179;

%Finding conditions:
if(z<7958)
    T1=Ts*(1-((gamma-1)/gamma)*(z/zstar));
    p1=ps*(1-((gamma-1)/gamma)*(z/zstar))^(gamma/(gamma-1));
else
    T1=210;
    p1=1000*33.6*exp((7959-z)/6605);
end
Tt1=T1*(1+m1^2*(gamma-1)/2);
pt1=p1*(1+m1^2*(gamma-1)/2)^(gamma/(gamma-1));
a1=sqrt(gamma*R*T1);
v1=m1*a1;
Cp1=a+b*T1;
%End of module 1
%%
%Module 2: For inlet/diffuser
gamma=1.4;
Tt2=Tt1;
T2=Tt2/(1+m2^2*(gamma-1)/2);
pt2=p1*(1+nd*m1^2*(gamma-1)/2)^(gamma/(gamma-1));
p2=pt2/(1+m2^2*(gamma-1)/2)^(gamma/(gamma-1));
Cp2=a+b*T2;
s12=Cp2*log(Tt2/Tt1)-R*log(pt2/pt1);
a2=sqrt(gamma*R*T2);
v2=m2*a2;
%End of module 2
%%
%Module 3: For combustor
gamma=1.3;
Tt3choked = Tt2*((1/(2*(gamma+1)))*(1/(m2^2))*((1+gamma*m2^2)^2)*(1+((gamma-1)/2)*m2^2)^-1);
%Checking thermal choking:
if(Tt3choked < Tt3max)
    m3=1;
    Tt3 = Tt3choked;
else
    Tt3 = Tt3max;
    C=(Tt3/Tt2)*(m2^2)*((1 + ((gamma - 1)/2)*m2^2)/((1 + gamma*(m2^2))^2));
    A=C*gamma^2-(gamma-1)/2;
    B=2*C*gamma-1;
    m3=sqrt((-1*B+sqrt(B^2-4*A*C))/(2*A));
    if(m3>m2 && m2<1 && m3<=1)
        
    else
        m3=sqrt((-1*B-sqrt(B^2-4*A*C))/(2*A));
    end
end
q23=a*(Tt3-Tt2)+0.5*b*(Tt3^2-Tt2^2);
p3=p2;
T3=Tt3/(1+m3^2*(gamma-1)/2);
pt3=p3*(1+m3^2*(gamma-1)/2)^(gamma/(gamma-1));
Cp3=a+b*T3;
a3=sqrt(gamma*R*T3);
v3=m3*a3;
s23=Cp3*log(Tt3/Tt2)-R*log(pt3/pt2);
%End of module 3
%%
%Module 4: Converging nozzle
gamma=1.3;
mtest=sqrt((2/(gamma-1))*(nn*(1-(p1/pt3)^((gamma-1)/gamma))/(1-nn*(1-(p1/pt3)^((gamma-1)/gamma)))));
if(mtest<1)
    me=mtest;
    pe=p1;
else
    me=1;
    pe=pt3*(1-nn^-1*((gamma-1)/(gamma+1)))^(gamma/(gamma-1));
end
Tte=Tt3;
Te=Tte/(1+me^2*(gamma-1)/2);
pte=pe*(1+me^2*(gamma-1)/2)^(gamma/(gamma-1));
ae=sqrt(gamma*R*Te);
ve=me*ae;
Cpe=a+b*Te;
exmassflux=pe*ve*Ae/(R*Te);
s3e=Cpe*log(Tte/Tt3)-R*log(pte/pt3);
se1=s3e+s23+s12;
%End of module 4
%%
%Module 5: External Flow
gamma=1.3;
%Checking nozzle choking condition:
if(mtest<1)
    nnext=1;
else
    nnext=mtest^(-0.3);
end
%Calculations
Tt4=Tte;
T4=Tt4*(1-nnext*(1-(p1/pte)^((gamma-1)/gamma)));
m4=sqrt((2/(gamma-1))*(Tt4/T4-1));
p4=p1;
pt4=p4*(1+m4^2*(gamma-1)/2)^(gamma/(gamma-1));
a4=sqrt(gamma*R*T4);
v4=m4*a4;
Cp4=a+b*T4;
s4e=Cp4*log(Tt4/Tte)-R*log(pt4/pte);
s41=se1+s4e;
%End of module 5
%%
%Module 6: Performance parameters
inmassflux=exmassflux/(1+q23/qf);
fuelmassflux=(exmassflux-inmassflux);
f=fuelmassflux/inmassflux;
jThrust=inmassflux*((1+f)*ve-v1);
pThrust=(pe-p1)*Ae;
totThrust=jThrust+pThrust;
TSFC=(fuelmassflux*3600)/totThrust;
Isp=totThrust/(fuelmassflux*9.8);
veq=ve+((pe-p1)*(Ae/exmassflux));
nth=((exmassflux*veq*veq/2)-(inmassflux*v1*v1/2))/(inmassflux*q23);
nprop=2/(1+(veq/v1));
noverall=nth*nprop;
propP=totThrust*v1;
output=totThrust;
end

%End of module 6
