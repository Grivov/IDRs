clear
clc

Folder='InitialState/';
mkdir(Folder)

Chain=450;
Spacing=8;
Sticker=ceil(Chain/2/Spacing)*2;
ParticleSize=3;

NA=Chain+Sticker; %all beads in a polymer
NPol=50; %number of polymers
NM=NPol*NA; %total polymer beads in a box
% NPar=10^4; %total particles in a box
% BoxSize(1)=400;

NPar=50; %total particles in a box
BoxSize(1)=100;
BoxSize(2)=25;
BoxSize(3)=25;

Replicates=30;
for rep=1:Replicates
    load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle' num2str(ParticleSize) 'nm.mat']);
    load(['InitialState/IdealChain_Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPol) '_Rep' num2str(rep) '.mat']);
    mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPol) '_Particle' num2str(NPar)];

    NB=3;
    x1=((1:(2*NB))-(2*NB+1)/2)*(BoxSize(2)/NB); %[-1 1]/2*(Template(end)-Template(1)+4);
    y1=((1:NB)-(NB+1)/2)*(BoxSize(2)/NB);
    z1=((1:NB)-(NB+1)/2)*(BoxSize(3)/NB);

    [X1,Y1,Z1]=meshgrid(x1,y1,z1);
    X1=reshape(X1,[],1);
    Y1=reshape(Y1,[],1);
    Z1=reshape(Z1,[],1);
    
    Index=randperm(length(X1));
    X1=X1(Index);
    Y1=Y1(Index);
    Z1=Z1(Index);

    Polymer1=zeros(3,NA,NPol);
    for np=1:NPol
        Polymer1(1,:,np)=XPol(np,:)+X1(np,1);
        Polymer1(2,:,np)=YPol(np,:)+Y1(np,1);
        Polymer1(3,:,np)=ZPol(np,:)+Z1(np,1);
    end

    atype=zeros(NA,1)+1;
    atype(Chain+(1:Sticker/2),1)=2;
    atype(Chain+Sticker/2+(1:Sticker/2),1)=3;

    figure(rep)
    for np=1:NPol
  	    plot3(Polymer1(1,atype==1,np),Polymer1(2,atype==1,np),Polymer1(3,atype==1,np),'k-'); hold on
        plot3(Polymer1(1,atype==2,np),Polymer1(2,atype==2,np),Polymer1(3,atype==2,np),'r.');
        plot3(Polymer1(1,atype==3,np),Polymer1(2,atype==3,np),Polymer1(3,atype==3,np),'b.'); 
    end
    axis equal
    axis([-BoxSize(1)/2 BoxSize(1)/2 -BoxSize(2)/2 BoxSize(2)/2 -BoxSize(3)/2 BoxSize(3)/2])

    if max(max(Polymer1(1,:,:))>29 | min(min(Polymer1(1,:,:))))<-29
        ['polymers outside the confinement box']
    end
    NB=Chain-1+Sticker;
    btype=zeros(NB,1);
    bond=zeros(NB,2);

    btype(1:(Chain-1),1)=1;
    btype(Chain:end,1)=2;
    bond(1:(Chain-1),1)=1:(Chain-1);
    bond(1:(Chain-1),2)=2:Chain;
    bond(Chain:end,:)=bondtype2';

    Monomer=zeros(3,NM);
    NBond=NB*NPol;
    Bond=zeros(2,NBond);
    Btype=zeros(1,NBond);
    Atype=zeros(1,NM);

    cut=30;
    x2=(-BoxSize(1)/2+BeadSize(end)/2):BeadSize(end):(BoxSize(1)/2-BeadSize(end)/2);
    x2=x2(x2>cut|x2<-cut);
    y2=(-BoxSize(2)/2+BeadSize(end)/2):BeadSize(end):(BoxSize(2)/2-BeadSize(end)/2);
    z2=(-BoxSize(3)/2+BeadSize(end)/2):BeadSize(end):(BoxSize(3)/2-BeadSize(end)/2);
    [X2,Y2,Z2]=meshgrid(x2,y2,z2);
    X2=reshape(X2,[],1);
    Y2=reshape(Y2,[],1);
    Z2=reshape(Z2,[],1);

    Particle=zeros(NPar,3);
    Index=randperm(length(X2));
    X2=X2(Index);
    Y2=Y2(Index);
    Z2=Z2(Index);
    for i=1:NPar
        Particle(i,1) = X2(i,1);
        Particle(i,2) = Y2(i,1);
        Particle(i,3) = Z2(i,1);
    end
    plot3(Particle(:,1),Particle(:,2),Particle(:,3),'go');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nm=0;
    nb=0;

    for np=1:NPol
        for n=1:NA
       	    nm=nm+1;
       	    Monomer(1,nm)=Polymer1(1,n,np);
       	    Monomer(2,nm)=Polymer1(2,n,np);
       	    Monomer(3,nm)=Polymer1(3,n,np);
            Atype(1,nm)=atype(n,1);
        end
        for n=1:NB
            nb=nb+1;
            Bond(1,nb)=bond(n,1)+(np-1)*NA;
            Bond(2,nb)=bond(n,2)+(np-1)*NA;
            Btype(1,nb)=btype(n,1);
        end
    end

    for i=1:NPar
        nm=nm+1;
        Monomer(1,nm)=Particle(i,1);
        Monomer(2,nm)=Particle(i,2);
        Monomer(3,nm)=Particle(i,3);
        Atype(1,nm)=4;
    end

    save([Folder mode '_Rep' num2str(rep) '.mat'],'Monomer','Bond','Atype','Btype','BoxSize','NA');
end
