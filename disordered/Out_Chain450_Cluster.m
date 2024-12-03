clear all
clc

Chain=450; 
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
Replicates=1;
for A=9
    ReadFolder=['StickerSpacer_Chain' num2str(Chain) '/Out/'];
    SaveFolder=['StickerSpacer_Chain' num2str(Chain) '/Out_ClusterAnalysis/'];
    mkdir(SaveFolder)

    for rep=1:Replicates
        ReadFilename=[ReadFolder mode '_A' num2str(A) '_Rep' num2str(rep) '.xyz'];
        [X,Y,Z]=ReadPosition(ReadFilename,NAtom,NT);
        for nt=1:NT
            IDden=1:NPolymerBeads;
            if std(X(nt,IDden))>30
                print('condensate split')
            end
	        while abs(mean(X(nt,IDden)))>0.1
                X(nt,:)=X(nt,:)-mean(X(nt,IDden));
      	        X(nt,:)=(X(nt,:)-BoxSize(1)).*(X(nt,:)>=BoxSize(1)/2)+...
               	        (X(nt,:)+BoxSize(1)).*(X(nt,:)<=-BoxSize(1)/2)+...
               	        X(nt,:).*(X(nt,:)>-BoxSize(1)/2).*(X(nt,:)<BoxSize(1)/2); 
            end
            [nt mean(X(nt,IDden))]
        end
        save([SaveFolder mode '_A' num2str(A) '_Rep' num2str(rep) '.mat'],'X','Y','Z');
    end
end

function [X,Y,Z]=ReadPosition(ReadFilename,NP,NT)
X=zeros(NT,NP);
Y=zeros(NT,NP);
Z=zeros(NT,NP);
fid=fopen(ReadFilename,'r');
NH=2; 
% nt=0;
% while ~feof(fid)
%     nt=nt+1;
for nt=1:NT
    nt
  	for nh=1:NH
      	Header=fgetl(fid);
    end
   	Data=textscan(fid,'%f %f %f %f',NP);
  	Header=fgetl(fid);
   	X(nt,:)=Data{1,2}; 
  	Y(nt,:)=Data{1,3};
  	Z(nt,:)=Data{1,4};
end
fclose(fid);
end
