clear all
%clc

Chain=450; 
NPol=50; 
NPar=1;
load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle3nm.mat']);
Spacing=8;
Sticker=ceil(Chain/2/Spacing)*2;
NPolymer=50; 
NParticle=1;
NPolymerBeads=NPolymer*(Sticker+Chain);
NAtom=NParticle+NPolymerBeads;
NT=50;

mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPolymer) '_Particle' num2str(NParticle)];
Folder='InitialState/';
load([Folder mode '.mat']);
load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle3nm.mat']);

NI=1000;
NK=100;
Cut=5;

PRadius=0;
Replicates=1;

for A=9

    ReadFolder=['StickerSpacer_Chain' num2str(Chain) '/Out_ClusterAnalysis/'];
    Probability1=zeros(NT,Replicates);
%     Probability2=zeros(NT,Replicates);
%     Probability3=zeros(NT,Replicates);
%     Probability123=zeros(NT,Replicates);
    VolumeFraction=zeros(Replicates,1);
    PartitionRatio=zeros(Replicates,1);

    for rep=1:Replicates
        load([ReadFolder mode '_A' num2str(A) '_Rep' num2str(rep) '.mat']);
        X1=X(:,Atype==1);
        Y1=Y(:,Atype==1);
        Z1=Z(:,Atype==1);
        [XE1,YE1,ZE1]=ExtendPosition(X1,Y1,Z1,BoxSize);
%         X2=X(:,Atype==2);
%         Y2=Y(:,Atype==2);
%         Z2=Z(:,Atype==2);
%         [XE2,YE2,ZE2]=ExtendPosition(X2,Y2,Z2,BoxSize);
%         X3=X(:,Atype==3);
%         Y3=Y(:,Atype==3);
%         Z3=Z(:,Atype==3);
%         [XE3,YE3,ZE3]=ExtendPosition(X3,Y3,Z3,BoxSize);
%         X4=X(:,Atype==4);
%         Y4=Y(:,Atype==4);
%         Z4=Z(:,Atype==4);
%         [XE4,YE4,ZE4]=ExtendPosition(X4,Y4,Z4,BoxSize);

        for nt=1:NT
            BX=40/2;
            BY=BoxSize(2)/2;
            BZ=BoxSize(3)/2;
            XI=(rand(NI,NT)-0.5)*2*BX;
  	        YI=(rand(NI,NT)-0.5)*2*BY;
  	        ZI=(rand(NI,NT)-0.5)*2*BZ;
  	        RI=[XI(:,nt),YI(:,nt),ZI(:,nt)]; 
            Keep1=XE1(:,nt)<(BX+Cut)&XE1(:,nt)>-(BX+Cut)&YE1(:,nt)<(BY+Cut)&YE1(:,nt)>-(BY+Cut)&ZE1(:,nt)<(BZ+Cut)&ZE1(:,nt)>-(BZ+Cut);
            XP1=XE1(Keep1,nt);
            YP1=YE1(Keep1,nt);
            ZP1=ZE1(Keep1,nt);
%             Keep2=XE2(:,nt)<(BX+Cut)&XE2(:,nt)>-(BX+Cut)&YE2(:,nt)<(BY+Cut)&YE2(:,nt)>-(BY+Cut)&ZE2(:,nt)<(BZ+Cut)&ZE2(:,nt)>-(BZ+Cut);
%             XP2=XE2(Keep2,nt);
%             YP2=YE2(Keep2,nt);
%             ZP2=ZE2(Keep2,nt);
%             Keep3=XE3(:,nt)<(BX+Cut)&XE3(:,nt)>-(BX+Cut)&YE3(:,nt)<(BY+Cut)&YE3(:,nt)>-(BY+Cut)&ZE3(:,nt)<(BZ+Cut)&ZE3(:,nt)>-(BZ+Cut);
%             XP3=XE3(Keep3,nt);
%             YP3=YE3(Keep3,nt);
%             ZP3=ZE3(Keep3,nt);
%             Keep4=XE4(:,nt)<(BX+Cut)&XE4(:,nt)>-(BX+Cut)&YE4(:,nt)<(BY+Cut)&YE4(:,nt)>-(BY+Cut)&ZE4(:,nt)<(BZ+Cut)&ZE4(:,nt)>-(BZ+Cut);
%             XP4=XE4(Keep4,nt);
%             YP4=YE4(Keep4,nt);
%             ZP4=ZE4(Keep4,nt);

  	        RP1=[XP1,YP1,ZP1];   
  	        [NList1,DList1]=knnsearch(RP1,RI,'K',NK);
            eps=1;
            sigma=BeadSize(1)/2+PRadius;
            if max((DList1(:,1)>(0.75*sigma)).*(DList1(:,NK)<(2^(1/6)*sigma)))==1
                'NK1 is not enough'
            end
            PE1=sum(LJ(DList1,eps,sigma),2);
            Probability1(nt,rep)=1-mean(exp(-PE1));

%   	      RP2=[XP2,YP2,ZP2];   
%   	      [NList2,DList2]=knnsearch(RP2,RI,'K',NK);
%             eps=1;
%             sigma=BeadSize(2)/2+PRadius;
%             if max((DList2(:,1)>(0.75*sigma)).*(DList2(:,NK)<(2^(1/6)*sigma)))==1
%                 'NK2 is not enough'
%             end
%             PE2=sum(LJ(DList2,eps,sigma),2);
%             Probability2(nt,rep)=1-mean(exp(-PE2));
% 
%   	      RP3=[XP3,YP3,ZP3];   
%   	      [NList3,DList3]=knnsearch(RP3,RI,'K',NK);
%             eps=1;
%             sigma=BeadSize(3)/2+PRadius;
%             if max((DList3(:,1)>(0.75*sigma)).*(DList3(:,NK)<(2^(1/6)*sigma)))==1
%                 'NK3 is not enough'
%             end
%             PE3=sum(LJ(DList3,eps,sigma),2);
%             Probability3(nt,rep)=1-mean(exp(-PE3));

%   	      RP4=[XP4,YP4,ZP4];   
%             NK4=min(NK,length(XP4));
%   	      [NList4,DList4]=knnsearch(RP4,RI,'K',NK4);
%             eps=1;
%             sigma=BeadSize(4)/2+PRadius;
%             if max((DList4(:,1)>(0.75*sigma)).*(DList4(:,NK4)<(2^(1/6)*sigma)))==1
%                 'NK4 is not enough'
%             end
%             PE4=sum(LJ(DList4,eps,sigma),2);
%             Probability4(nt,rep)=1-mean(exp(-PE4));

        end
        VolumeFraction(rep,1)=mean(mean(Probability1));
        PartitionRatio(rep,1)=1/(1-mean(mean(Probability1)));
    end
end

['BeadRadius ' num2str(PRadius) ' VolumeFraction ' num2str(mean(VolumeFraction))]
['BeadRadius ' num2str(PRadius) ' Partition ratio ' num2str(mean(PartitionRatio))]


function PE=LJ(r,eps,sigma)
PE=4*eps*((sigma./r).^12-(sigma./r).^6+1/4).*(r<=sigma*2^(1/6));
end

function [XE,YE,ZE]=ExtendPosition(X,Y,Z,BoxSize)
[NT,NP]=size(X);
[XS,YS,ZS]=meshgrid(0,[-1 0 1]*BoxSize(2),[-1 0 1]*BoxSize(3));
XS=reshape(XS,[],1);
YS=reshape(YS,[],1);
ZS=reshape(ZS,[],1);
XE=zeros(9*NP,NT);
YE=zeros(9*NP,NT);
ZE=zeros(9*NP,NT);
for n=1:9
   	XE((n-1)*NP+(1:NP),:)=X'+XS(n);
   	YE((n-1)*NP+(1:NP),:)=Y'+YS(n);
  	ZE((n-1)*NP+(1:NP),:)=Z'+ZS(n);
end
end