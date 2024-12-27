clear all
%clc

R=0.6/2;
H=R-0.38/2;
V1=pi/3*(3*R-H)*H^2;
V2=4*pi/3*R^3-2*V1;

Chain=10; 
Sticker=10;
NPolymer=400; 
NParticle=1;
NPolymerBeads=NPolymer*(Sticker+Chain);
NAtom=NParticle+NPolymerBeads;
NT=50;

mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPolymer) '_Particle' num2str(NParticle)];
Folder='InitialState/';
load([Folder mode '.mat']);

load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle3nm.mat'])
DX=BoxSize(1)/50;
XSize=-(BoxSize(1)/2-DX/2):DX:(BoxSize(1)/2-DX/2);

Replicates=10;
for A=9.5
    ReadFolder=['StickerSpacer_Chain' num2str(Chain) '/Out_ClusterAnalysis/'];
    for rep=6%1:Replicates
        load([ReadFolder mode '_A' num2str(A) '_Rep' num2str(rep) '.mat']);
        XPol=reshape(X(:,1:NPolymerBeads),[],1);
        XPar=reshape(X(:,(NPolymerBeads+1):end),[],1);
        CountPol=hist(XPol,XSize)/NT;
        CountPar=hist(XPar,XSize)/NT;
        ConcentrationPol=CountPol/(6.02*10^23)/(BoxSize(2)*BoxSize(3)*DX*10^-27); %mM
        ConcentrationPar=CountPar/(6.02*10^23)/(BoxSize(2)*BoxSize(3)*DX*10^-27); %mM
        VolumeFracPol=4*pi/3*(BeadSize(1)/2)^3*Chain/(Sticker+Chain)*CountPol/(BoxSize(2)*BoxSize(3)*DX); 
        VolumeFracPar=CountPar*4*pi/3*(BeadSize(4)/2)^3/(BoxSize(2)*BoxSize(3)*DX); 
    end

figure(1)
semilogy(XSize,ConcentrationPol/(Sticker+Chain),'s-'); hold on
%plot(XSize,ConcentrationPar,'o-'); hold on
figure(2)
plot(XSize,VolumeFracPol,'s-'); hold on
plot(XSize,VolumeFracPar,'s-'); 
plot(XSize,VolumeFracPar+VolumeFracPol,'s-'); 
xlabel('x (nm)')
ylabel('Volume Fraction')
legend('Polymer','Particle','Total')
legend boxoff

end
% figure(3)
% s=1;
% logC=log(ConcentrationPol);
% plot(XSize,(-logC+max(logC))/s,'o-'); hold on
% plot(XSize,0.0025*(XSize-63).^2,'k-');
% plot(XSize,0.0025*(XSize+63).^2,'k-');
% axis([-150 150 0 20])