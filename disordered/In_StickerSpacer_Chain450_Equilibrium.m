clear
clc

% Protocol 2*10^7  anneal soft potential
%          2*10^7  remove wall
%          16*10^7 equilibration
%          10^8    record data 400 steps
%          5*10^7  remove particle equilibration
%          10^8    record data 400 steps

global Monomer Bond Atype Btype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp RelaxationTime NA

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

Chain=450; 
NPol=50; 
NPar=50;
load(['Parameter/Parameter_Chain' num2str(Chain) '_Particle3nm.mat']);
Spacing=8;
Sticker=ceil(Chain/2/Spacing)*2;
mode=['Sticker' num2str(Sticker) '_Chain' num2str(Chain) '_NP' num2str(NPol) '_Particle' num2str(NPar)];
Folder='InitialState/';

InFolder=['StickerSpacer_Chain' num2str(Chain) '/'];
mkdir(InFolder)
OutFolder=['Out/'];
mkdir([InFolder OutFolder]);

index=0;
Replicates=30;
for A=9
    for rep=1:Replicates
        index=index+1;
        load([Folder mode '_Rep' num2str(rep) '.mat']);
        InitCondFilename=[mode '_Rep' num2str(rep) '.initial'];
        InitCondGenerate(InFolder,InitCondFilename);
        InFilename=['Equilibrium_' mode '_Index_' num2str(index) '.in'];
        OutFilename=[OutFolder mode '_A' num2str(A) '_Rep' num2str(rep)];
        InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,A)
    end
end

%%
function []=InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,A)
global BeadSize Temp kBT Damp

eps=kBT;
sigma1=BeadSize(1);
sigma2=BeadSize(2);
sigma3=BeadSize(3);
sigma4=BeadSize(4);
sigma11=(sigma1+sigma1)/2;
sigma12=(sigma1+sigma2)/2;
sigma13=(sigma1+sigma3)/2;
sigma14=(sigma1+sigma4)/2;
sigma22=(sigma2+sigma2)/2;
sigma23=(sigma2+sigma3)/2;
sigma24=(sigma2+sigma4)/2;
sigma33=(sigma3+sigma3)/2;
sigma34=(sigma3+sigma4)/2;
sigma44=(sigma4+sigma4)/2;

A=-A*kBT;
bondlength1=0.38;
bondlength2=0.3;
Rb=30;

TimeStep=Damp/100;
RunSteps=10^7;
Thermo=RunSteps;
RecordSteps=RunSteps/10;

fid=fopen([InFolder InFilename],'w');
fprintf(fid, ['units nano\n']);
fprintf(fid, ['boundary p p p\n']);
fprintf(fid, ['atom_style bond\n\n']);

%fprintf(fid, ['processors 2 * *\n\n']);

fprintf(fid, ['read_data ' InitCondFilename '\n\n']);

s=3;
fprintf(fid, ['pair_style hybrid lj/cut ' num2str(sigma44*2^(1/6)) ' soft ' num2str(sigma11*2^(1/6)) ' zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff * * zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff 1 1 soft ' num2str(s*eps) ' ' num2str(sigma11*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 2 soft ' num2str(s*eps) ' ' num2str(sigma22*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 3 3 soft ' num2str(s*eps) ' ' num2str(sigma33*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 3 soft ' num2str(0) ' ' num2str(sigma23/2) '\n']);
fprintf(fid, ['pair_coeff 1 4 lj/cut ' num2str(eps) ' ' num2str(sigma14) ' ' num2str(sigma14*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 4 4 lj/cut ' num2str(eps) ' ' num2str(sigma44) ' ' num2str(sigma44*2^(1/6)) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 0.0 1.0 1.0\n\n']);

fprintf(fid, ['bond_style harmonic\n']);
fprintf(fid, ['bond_coeff 1 ' num2str(100*kBT/bondlength1^2) ' ' num2str(bondlength1) '\n']);
fprintf(fid, ['bond_coeff 2 ' num2str(100*kBT/bondlength2^2) ' ' num2str(bondlength2) '\n\n']);

%fprintf(fid, ['neighbor ' num2str(5-sigma44*2^(1/6)) ' bin\n']);

fprintf(fid, ['neighbor ' num2str(5-sigma44*2^(1/6)) ' multi\n\n']);
fprintf(fid, ['comm_style tiled\n']);
fprintf(fid, ['neigh_modify every 1 delay 0\n\n']);

fprintf(fid, ['region wallx block -' num2str(Rb) ' ' num2str(Rb) ...
                                ' -' num2str(30) ' ' num2str(30) ...
                                ' -' num2str(30) ' ' num2str(30) ' open 3 open 4 open 5 open 6\n\n']);

fprintf(fid, ['group polymer type 1 2 3\n']);
fprintf(fid, ['group particle type 4\n\n']);

fprintf(fid, ['fix 1 all nve\n']);
fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero yes\n']);
fprintf(fid, ['fix 4 polymer wall/region wallx lj126 ' num2str(eps) ' ' num2str(sigma11) ' ' num2str(sigma11*2^(1/6)) '\n']);
%fprintf(fid, ['fix 5 all balance ' num2str(RunSteps/100) ' 1.05 shift x 10 1.05\n\n']);
fprintf(fid, ['fix 5 all balance ' num2str(RunSteps/100) ' 1.1 rcb\n\n']);

fprintf(fid, ['thermo ' num2str(Thermo) '\n']);
fprintf(fid, ['timestep ' num2str(TimeStep) '\n\n']);

fprintf(fid, ['dump 1 all movie ' num2str(RecordSteps) ' ' OutFilename '.mpeg type type zoom 4 box yes 0.01 bond type 0.2 view 88 88 size 1000 400 shiny 0.5\n']);
fprintf(fid, ['dump_modify 1 acolor 1 white\n']); 
fprintf(fid, ['dump_modify 1 acolor 2 red\n']);
fprintf(fid, ['dump_modify 1 acolor 3 blue\n']);
fprintf(fid, ['dump_modify 1 acolor 4 green\n']);
fprintf(fid, ['dump_modify 1 bcolor 1 white\n']);
fprintf(fid, ['dump_modify 1 bcolor 2 white\n']);
fprintf(fid, ['dump_modify 1 adiam 1 ' num2str(sigma11) '\n']); 
fprintf(fid, ['dump_modify 1 adiam 2 ' num2str(sigma22) '\n']); 
fprintf(fid, ['dump_modify 1 adiam 3 ' num2str(sigma33) '\n']);
fprintf(fid, ['dump_modify 1 adiam 4 ' num2str(sigma44) '\n\n']);

fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);

s=50;
fprintf(fid, ['pair_style hybrid lj/cut ' num2str(sigma44*2^(1/6)) ' soft ' num2str(sigma11*2^(1/6)) ' zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff * * zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff 1 1 soft ' num2str(s*eps) ' ' num2str(sigma11*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 2 soft ' num2str(s*eps) ' ' num2str(sigma22*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 3 3 soft ' num2str(s*eps) ' ' num2str(sigma33*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 3 soft ' num2str(0) ' ' num2str(sigma23/2) '\n']);
fprintf(fid, ['pair_coeff 1 4 lj/cut ' num2str(eps) ' ' num2str(sigma14) ' ' num2str(sigma14*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 4 4 lj/cut ' num2str(eps) ' ' num2str(sigma44) ' ' num2str(sigma44*2^(1/6)) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 0.0 1.0 1.0\n\n']);

fprintf(fid, ['run ' num2str(RunSteps/10) '\n\n']);

fprintf(fid, ['pair_style hybrid lj/cut ' num2str(sigma44*2^(1/6)) ' soft ' num2str(sigma11*2^(1/6)) ' zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff * * zero ' num2str(sigma11) '\n']);
fprintf(fid, ['pair_coeff 1 1 lj/cut ' num2str(eps) ' ' num2str(sigma11) ' ' num2str(sigma11*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 2 lj/cut ' num2str(eps) ' ' num2str(sigma22) ' ' num2str(sigma22*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 3 3 lj/cut ' num2str(eps) ' ' num2str(sigma33) ' ' num2str(sigma33*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 2 3 soft ' num2str(0) ' ' num2str(sigma23/2) '\n']);
fprintf(fid, ['pair_coeff 1 4 lj/cut ' num2str(eps) ' ' num2str(sigma14) ' ' num2str(sigma14*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 4 4 lj/cut ' num2str(eps) ' ' num2str(sigma44) ' ' num2str(sigma44*2^(1/6)) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 0.0 1.0 1.0\n\n']);

fprintf(fid, ['variable A equal "ramp(' num2str(0) ',' num2str(A) ')"\n']);
fprintf(fid, ['fix 3 all adapt 1 pair soft a 2 3 v_A\n\n']);
fprintf(fid, ['run ' num2str(RunSteps*4) '\n\n']);

fprintf(fid, ['unfix 3\n']);
fprintf(fid, ['pair_coeff 2 3 soft ' num2str(A) ' ' num2str(sigma11/2) '\n']);
fprintf(fid, ['run ' num2str(RunSteps*1) '\n\n']);

fprintf(fid, ['unfix 4\n']);
fprintf(fid, ['run ' num2str(RunSteps*5) '\n\n']);

fprintf(fid, ['dump 4 all xyz ' num2str(RecordSteps) ' ' OutFilename '.xyz\n']); 
fprintf(fid, ['run ' num2str(RunSteps*5) '\n\n']);

fprintf(fid, ['write_restart ' OutFilename '.restart\n\n']);

fclose(fid); 

end

function []=InitCondGenerate(InFolder,InitCondFilename)

global Monomer Bond Atype Btype BoxSize BeadMass NA

Natom_type=4; %number of atom types;
Natom=length(Monomer);
Nbond_type=length(unique(Btype)); %number of bond types;
Nbond=length(Bond); %number of DNA bonds;
Nangl_type=0; %number of angle types;
Nangl=0; %number of angles;
Ndihe_type=0; %number of dihedral types;
Ndihe=0; %number of dihedrals;
Nimpr_type=0; %number of improper types;
Nimpr=0; %number of impropers;

xlo=-BoxSize(1)/2; xhi=BoxSize(1)/2; %x boundary
ylo=-BoxSize(2)/2; yhi=BoxSize(2)/2; %y boundary
zlo=-BoxSize(3)/2; zhi=BoxSize(3)/2; %z boundary

V=zeros(3,Natom);

fid=fopen([InFolder InitCondFilename],'w');
fprintf(fid,'LAMMPS chain data file\n\n');
fprintf(fid,'%d atoms\n', Natom);
fprintf(fid,'%d bonds\n', Nbond);
fprintf(fid,'%d angles\n', Nangl);
fprintf(fid,'%d dihedrals\n', Ndihe);
fprintf(fid,'%d impropers\n\n', Nimpr);
fprintf(fid,'%d atom types\n', Natom_type);
fprintf(fid,'%d bond types\n', Nbond_type);
fprintf(fid,'%d angle types\n', Nangl_type);
fprintf(fid,'%d dihedral types\n', Ndihe_type);
fprintf(fid,'%d improper types\n\n', Nimpr_type);
fprintf(fid,'%8.5f %8.5f xlo xhi\n', xlo, xhi);
fprintf(fid,'%8.5f %8.5f ylo yhi\n', ylo, yhi);
fprintf(fid,'%8.5f %8.5f zlo zhi\n\n', zlo, zhi);

fprintf(fid,'Masses\n\n');
for n=1:Natom_type
    fprintf(fid,'%d %8.5f\n',n,BeadMass(n));
end
fprintf(fid,'\n');

Mole_type=1;
fprintf(fid,'Atoms\n\n');
for i=1:Natom
    Atom_type=Atype(i);
    fprintf(fid,[num2str(i) ' ' num2str(Mole_type) ' ' num2str(Atom_type) ' ' ...
                 num2str(Monomer(1,i)) ' ' num2str(Monomer(2,i)) ' ' num2str(Monomer(3,i)) '\n']);
    if mod(i,NA)==0
        Mole_type=Mole_type+1;
    end    
end
fprintf(fid,'\n');

fprintf(fid,'Velocities\n\n');
for i=1:Natom
    fprintf(fid,[num2str(i) ' ' num2str(V(1,i)) ' ' num2str(V(2,i)) ' ' num2str(V(3,i)) '\n']);
end
fprintf(fid,'\n');

fprintf(fid,'Bonds\n\n');
for i=1:Nbond
    fprintf(fid,'%d %d %d %d\n',...
            i,Btype(i),Bond(1,i),Bond(2,i));
end
fprintf(fid,'\n');

fclose(fid); 
end



function DeltaF=A2DeltaF(A,rcut)
dr=10^-4;
r=dr/2:dr:(rcut-dr/2);
DeltaF=zeros(length(A),1);
for na=1:length(A)
    DeltaF(na,1)=log(sum(4*pi*r.^2.*exp(U_Soft(r,A(na),rcut)))*dr/sum(4*pi*r.^2*dr));
end
end

function y=U_Soft(r,A,rcut)
y=A*(1+cos(r*pi/rcut)).*(r<rcut);
end