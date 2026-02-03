# Isothermal-reactor-model
%%%% Isothemal Fixed Bed Reactor
function dF=isothermal(w,F)
dF=zeros(9,1);
%%%Parameters
R=8.314;                              % ideal gas constant
T=806;                                % Temperature (K)
p=1.25;                               % total pressure (bar)
epslon=0.4312;                        % void fraction
rou=1422;                             % bulk density (kg/m^3)
r=3.5;                                % radius of reactor (m)
Area=pi*r^2;                          % area (m^2)
a=epslon/rou;                         % connversion factor (unoccupied volume m^3/mass Kg)
Keq=0.336;                            % Equilibrium constant (1/bar)
% Rate constant
k1t=2.2215*10^16*exp(-272230/(R*T));  % rate constant first thermal reaction (kmol/m^3.hr.bar)
k2t=2.4217*10^20*exp(-352790/(R*T));  % rate constant second thermal reaction (kmol/m^3.hr.bar)
k3t=3.8224*10^17*exp(-313060/(R*T));  % rate constant third thermal reaction (kmol/m^3.hr.bar)
k1=4.594*10^9*exp(-175380/(R*T));     % rate constant first catalytic reaction (kmol/(kgcat·hr)
k2=1.06*10^15*exp(-296290/(R*T));     % rate constant second catalytic reaction (kmol/(kgcat·hr)
k3=1.246*10^26*exp(-474760/(R*T));    % rate constant third catalytic reaction (kmol/(kgcat·hr)
k4=8.2024*10^10*exp(-213780/(R*T));   % rate constant fourth catalytic reaction (kmol/(kgcat·hr)
KEB=1.014*10^-5*exp(102220/(R*T));    % adsorption constant of ethylbenzene (1/bar)
KST=2.678*10^-5*exp(104560/(R*T));    % adsorption constant of Styrene (1/bar)
KH=4.519*10^-7*exp(117950/(R*T));     % adsorption constant of hydrogene (1/bar)
%flow rate in terms of pressure
Ft=F(1)+F(2)+F(3)+F(4)+F(5)+F(6)+F(7)+F(8)+F(9); %total flow rate (kmol/hr)
yEB=F(1)/Ft;                          % mole fraction of ethylbenzene
PEB=yEB*p;                            % partial pressure of ethylbenzene (bar)
yST=F(2)/Ft;                          % mole fraction of Styrene
pST=yST*p;                            % partial pressure of Styrene (bar)
yBZ=F(3)/Ft;                          % mole fraction of benzene
pBZ=yBZ*p;                            % partial pressure of benzene (bar)
yTo=F(4)/Ft;                          % mole fraction of toluene
pTo=yTo*p;                            % partial pressure of toluene (bar)
yH=F(5)/Ft;                           % mole fraction of hydrogene
pH=yH*p;                              % partial pressure of hydrogene (bar)
yMeth=F(6)/Ft;                        % mole fraction of methane
pMeth=yMeth*p;                        % partial pressure of methane (bar)
yEth=F(7)/Ft;                         % mole fraction of ethylene
pEth=yEth*p;                          % partial pressure of ethylene (bar)
yN=F(8)/Ft;                           % mole fraction of Nitrogene
pN=yN*p;                              % partial pressure of Nitrogene (bar)
yH2O=F(9)/Ft;                         % mole fraction of steam
pH2O=yH2O*p;                          % partial pressure pf steam (bar)
%rate laws
rt1=k1t*(PEB-pST*pH/Keq);             % first thermal reaction (kmol/m^3.hr)
rt2=k2t*PEB;                          % second thermal reaction (kmol/m^3.hr)
rt3=k3t*PEB;                          % third thermal reaction (kmol/m^3.hr)
den=(1+KEB*PEB+KH*pH+KST*pST)^2;      % denominator adsorption term 
r1=k1*KEB*(PEB-pST*pH/Keq)/den;       % first catalytic reaction (kmol/(kgcat·hr)
r2=k2*KEB*PEB/den;                    % second catalytic reaction (kmol/(kgcat·hr)
r3=k3*KEB*PEB*KH*pH/den;              % Third catalytic reaction (kmol/(kgcat·hr)
r4=k4*KST*pST*KH*pH/den;              % fourth catalytic reaction (kmol/(kgcat·hr)
% diffrential Equations with net rate
dF(1,1)=-r1-r2-r3-rt1*a-rt2*a-rt3*a;  % Etylbenzene mole balance (kmol/(kgcat·hr)
dF(2,1)=r1-r4+rt1*a;                  % styrene mole balance (kmol/(kgcat·hr)
dF(3,1)=rt2*a+r2;                     % benzene mole balance (kmol/(kgcat·hr)
dF(4,1)=r3+r4+rt3*a;                  % Toluene mole balance (kmol/(kgcat·hr)
dF(5,1)=rt1*a-rt3*a+r1-r3-2*r4;       % hydrogene mole balance (kmol/(kgcat·hr)
dF(6,1)=rt3*a+r3+r4;                  % methane mole balance (kmol/(kgcat·hr)
dF(7,1)=r2+rt2*a;                     % Ethylene mole balance (kmol/(kgcat·hr)
dF(8,1)=0;                            % Nitrogene mole balance (kmol/(kgcat·hr)
dF(9,1)=0;                            % steam mole balance (kmol/(kgcat·hr)
%[w,F]=ode45(@isothermal,[0 72950],[707 7.104 0.293 4.968 0 0 0 0 7777])
%plotyy(w,F(:,1),w,F(:,2))
%legend('EB','ST')
%xlabel('weight of catalyst')
%ylabel('molar flow')
%conv=(707-F(:,1))/707*100;
%plot(w,conv)
%plot(w,conv)
%length=w(:,1)/(pi*3.5^2*1422);
%plot(length,conv)
end