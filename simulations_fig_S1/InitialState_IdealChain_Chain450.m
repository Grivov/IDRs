clear
clc

Folder='InitialState/';
mkdir(Folder)

Chain=450;
Spacing=8;
Sticker=ceil(Chain/2/Spacing)*2;
ParticleSize=4;
load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle' num2str(ParticleSize) 'nm.mat']);

NA=Sticker+Chain; %all beads in a polymer
NPol=50; %number of polymers
NM=NPol*NA; %total polymer beads in a box
% NPar=10^4; %total particles in a box
% BoxSize(1)=400;

mode=['IdealChain_Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPol)];

PolymerBond=0.38; %bond length
StickerBond=0.3;

Replicates=30;
for rep=1:Replicates
    % random walk
    XPol=zeros(NPol,NA);
    YPol=zeros(NPol,NA);
    ZPol=zeros(NPol,NA);
    for npol=1:NPol
        for n=2:Chain
            XPol(npol,n)=XPol(npol,n-1);
            YPol(npol,n)=YPol(npol,n-1);
            ZPol(npol,n)=ZPol(npol,n-1);
            s=randi(3);
            if s==1
                XPol(npol,n)=XPol(npol,n)+PolymerBond*(2*randi(2)-3);
            elseif s==2
                YPol(npol,n)=YPol(npol,n)+PolymerBond*(2*randi(2)-3);
            elseif s==3
                ZPol(npol,n)=ZPol(npol,n)+PolymerBond*(2*randi(2)-3);
            end
        end
    end

    Index1=1:Spacing:(Chain/2);
    Index2=Chain:(-Spacing):(Chain/2+1);
    Index2=Index2(end:-1:1);
    Sticker1X=XPol(:,Index1)+StickerBond;
    Sticker1Y=YPol(:,Index1);
    Sticker1Z=ZPol(:,Index1);
    Sticker2X=XPol(:,Index2)+StickerBond;
    Sticker2Y=YPol(:,Index2);
    Sticker2Z=ZPol(:,Index2);
    XPol(:,Chain+(1:Sticker/2))=Sticker1X;
    YPol(:,Chain+(1:Sticker/2))=Sticker1Y;
    ZPol(:,Chain+(1:Sticker/2))=Sticker1Z;
    XPol(:,Chain+Sticker/2+(1:Sticker/2))=Sticker2X;
    YPol(:,Chain+Sticker/2+(1:Sticker/2))=Sticker2Y;
    ZPol(:,Chain+Sticker/2+(1:Sticker/2))=Sticker2Z;
    XPol=XPol-mean(XPol,2);
    YPol=YPol-mean(YPol,2);
    ZPol=ZPol-mean(ZPol,2);
    
    bondtype2=[Index1,Index2;Chain+(1:Sticker/2),Chain+Sticker/2+(1:Sticker/2)];
    save([Folder mode '_Rep' num2str(rep) '.mat'],'XPol','YPol','ZPol','bondtype2');
end

npol=10;
plot3(XPol(npol,1:Chain),YPol(npol,1:Chain),ZPol(npol,1:Chain),'.-'); hold on
plot3(XPol(npol,Chain+(1:Sticker/2)),YPol(npol,Chain+(1:Sticker/2)),ZPol(npol,Chain+(1:Sticker/2)),'*');
plot3(XPol(npol,Chain+Sticker/2+(1:Sticker/2)),YPol(npol,Chain+Sticker/2+(1:Sticker/2)),ZPol(npol,Chain+Sticker/2+(1:Sticker/2)),'*');
axis equal
