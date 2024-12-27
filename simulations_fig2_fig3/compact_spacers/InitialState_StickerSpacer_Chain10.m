clear
clc

Folder='InitialState/';
mkdir(Folder)

Chain=10;
Side=1;
Sticker=Chain*Side;
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

load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle' num2str(ParticleSize) 'nm.mat']);
mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPol) '_Particle' num2str(NPar)];

Template=zeros(NA,3);
Template(1:Chain,1)=((1:Chain)-(Chain+1)/2)*BeadSize(1)*2;
Template(1*Chain+(1:Chain),1)=Template(1:Chain,1)+BeadSize(1)/2;
% Template(2*Chain+(1:Chain),1)=Template(1:Chain,1)-BeadSize(1)/2;
% Template(3*Chain+(1:Chain),1)=Template(1:Chain,1);
% Template(3*Chain+(1:Chain),2)=Template(1:Chain,2)+BeadSize(1)/2;
% Template(4*Chain+(1:Chain),1)=Template(1:Chain,1);
% Template(4*Chain+(1:Chain),2)=Template(1:Chain,2)-BeadSize(1)/2;
% Template(5*Chain+(1:Chain),1)=Template(1:Chain,1);
% Template(5*Chain+(1:Chain),3)=Template(1:Chain,3)+BeadSize(1)/2;
% Template(6*Chain+(1:Chain),1)=Template(1:Chain,1);
% Template(6*Chain+(1:Chain),3)=Template(1:Chain,3)-BeadSize(1)/2;

NB=10;
x1=0;
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
    Polymer1(1,:,np)=Template(:,1)+X1(np,1);
    Polymer1(2,:,np)=Template(:,2)+Y1(np,1);
    Polymer1(3,:,np)=Template(:,3)+Z1(np,1);
end

atype=zeros(NA,1)+3;
atype(1:Chain,1)=1;
for n=1:(Chain/2)
    atype((Chain+n):Chain:((Side+1)*Chain),1)=2;
end

mtype=zeros(NA,1);
for m=1:Chain
    mtype(m:Chain:((Side+1)*Chain),1)=m;
end

figure(1)
for np=1:NPol
  	plot3(Polymer1(1,atype==1,np),Polymer1(2,atype==1,np),Polymer1(3,atype==1,np),'ko-'); hold on
    plot3(Polymer1(1,atype==2,np),Polymer1(2,atype==2,np),Polymer1(3,atype==2,np),'r.');
    plot3(Polymer1(1,atype==3,np),Polymer1(2,atype==3,np),Polymer1(3,atype==3,np),'b.'); 
end
axis equal
axis([-BoxSize(1)/2 BoxSize(1)/2 -BoxSize(2)/2 BoxSize(2)/2 -BoxSize(3)/2 BoxSize(3)/2])

%%
NB=Chain-1+Sticker;
btype=zeros(NB,1);
bond=zeros(NB,2);

btype(1:(Chain-1),1)=1;
btype(Chain:end,1)=2;
bond(1:(Chain-1),1)=1:(Chain-1);
bond(1:(Chain-1),2)=2:Chain;
bond(1*Chain-1+(1:Chain),1)=1:Chain;
bond(1*Chain-1+(1:Chain),2)=1*Chain+(1:Chain);

Monomer=zeros(3,NM);
NBond=NB*NPol;
Bond=zeros(2,NBond);
Btype=zeros(1,NBond);
Atype=zeros(1,NM);
Mtype=zeros(1,NM);

cut=25;
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
        Mtype(1,nm)=mtype(n,1)+(np-1)*Chain;
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
    Mtype(1,nm)=0;
end

save([Folder mode '.mat'],'Monomer','Bond','Atype','Btype','Mtype','BoxSize','NA');
