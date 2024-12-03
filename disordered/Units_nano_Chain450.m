clear
clc

% mass = attograms = 10^-18 gram
% distance = nanometers
% time = nanoseconds
% energy = attogram-nanometer^2/nanosecond^2
% velocity = nanometers/nanosecond
% force = attogram-nanometer/nanosecond^2
% torque = attogram-nanometer^2/nanosecond^2
% temperature = Kelvin
% pressure = attogram/(nanometer-nanosecond^2)
% dynamic viscosity = attogram/(nanometer-nanosecond)
% charge = multiple of electron charge (1.0 is a proton)
% dipole = charge-nanometer
% electric field = volt/nanometer
% density = attograms/nanometer^dim

SaveFolder='Parameter/';
mkdir(SaveFolder);

Chain=450;
BeadSize1=0.6;
BeadSize2=0.6;
BeadSize3=0.6;
BeadSize4=3;
BeadSize=[BeadSize1,BeadSize2,BeadSize3,BeadSize4];

WaterEta=1; %ag/nm/ns %0.001 kg/m/s
BeadCsi=6*pi*WaterEta*BeadSize/2; %ag/ns

Temp=300; %K
kB=1.38*10^-2; %ag*nm^2/ns^2/K 1.38*10^-23 Kg*m^2/s^2/K
kBT=kB*Temp;
    
Damp=1;
D=kBT./BeadCsi;
BeadMass=Damp*BeadCsi;

dt=Damp/100;
T=dt:1:2*10^6;
loglog(T,(3*kBT/BeadMass(1)*T.^2).^0.5,'-'); hold on %1/2 m v^2 =3/2 kB T v^2*T^2
loglog(T,(3*kBT/BeadMass(2)*T.^2).^0.5,'-');
loglog(T,(3*kBT/BeadMass(3)*T.^2).^0.5,'-');
loglog(T,(6*D(1)*T).^0.5,'-');
loglog(T,(6*D(2)*T).^0.5,'-');
loglog(T,(6*D(3)*T).^0.5,'-');

save([SaveFolder 'Parameter_Chain' num2str(Chain) '_Particle' num2str(BeadSize(4)) 'nm.mat'],'BeadSize','BeadCsi','Damp','Temp','kBT','BeadMass');

