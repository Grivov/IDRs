clear all
%clc

R=0.6/2;
H=R-0.38/2;
V1=pi/3*(3*R-H)*H^2;
V2=4*pi/3*R^3-2*V1;

Chain=450; 
NPol=50; 
NPar=50;
load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle3nm.mat']);
Spacing=8;
Sticker=ceil(Chain/2/Spacing)*2;
NPolymer=50; 
NParticle=50;
NPolymerBeads=NPolymer*(Sticker+Chain);
NAtom=NParticle+NPolymerBeads;
NT=50;

mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPolymer) '_Particle' num2str(NParticle)];
Folder='InitialState/';
load(['InitialState/Sticker58_Chain450_NP50_Particle50_Rep1.mat']);

BeadSize=[0.6,0.6,0.6,3];
DX=BoxSize(1)/50;
XSize=-(BoxSize(1)/2-DX/2):DX:(BoxSize(1)/2-DX/2);
NX=length(XSize);

Replicates=1;
ConcentrationPol=zeros(Replicates,NX);
ConcentrationPar=zeros(Replicates,NX);
VolumeFracPol=zeros(Replicates,NX);
VolumeFracPar=zeros(Replicates,NX);
for A=9
    ReadFolder=['StickerSpacer_Chain' num2str(Chain) '/Out_ClusterAnalysis/'];
    for rep=1%1:Replicates
        load([ReadFolder mode '_A' num2str(A) '_Rep' num2str(rep) '.mat']);
        XPol=reshape(X(:,1:NPolymerBeads),[],1);
        XPar=reshape(X(:,(NPolymerBeads+1):end),[],1);
        CountPol=hist(XPol,XSize)/NT;
        CountPar=hist(XPar,XSize)/NT;
        ConcentrationPol(rep,:)=CountPol/(6.02*10^23)/(BoxSize(2)*BoxSize(3)*DX*10^-27); %mM
        ConcentrationPar(rep,:)=CountPar/(6.02*10^23)/(BoxSize(2)*BoxSize(3)*DX*10^-27); %mM
        VolumeFracPol(rep,:)=V2*Chain/(Sticker+Chain)*CountPol/(BoxSize(2)*BoxSize(3)*DX); 
        VolumeFracPar(rep,:)=CountPar*4*pi/3*(BeadSize(4)/2)^3/(BoxSize(2)*BoxSize(3)*DX); 
    end
end
figure(1)
semilogy(XSize,mean(ConcentrationPol,1)/(Sticker+Chain),'s-'); hold on
%plot(XSize,mean(ConcentrationPar,1),'o-'); hold on
xlabel('x (nm)')
ylabel('Polymer concentration (mM)')

figure(2)
plot(XSize,mean(VolumeFracPol,1),'s-'); hold on
plot(XSize,mean(VolumeFracPar,1),'s-'); 
plot(XSize,mean(VolumeFracPar+VolumeFracPol,1),'s-'); 
xlabel('x (nm)')
ylabel('Volume Fraction')
legend('Polymer','Particle','Total')
legend boxoff

mean(VolumeFracPol(:,XSize<20&XSize>-20))
