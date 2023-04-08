clear all
close all

% types
Natomtypes=6;
Nbondtypes=3;
Nangletypes=3;
Ndihedraltypes=3;  
Nimpropertypes=2;
cptmol=0;
cptatom=0; 
cptbond=0; 
cptangle=0; 
cptdihedrals=0; 
cptimpropers=0; 

% FeO walls
dFeO=2.085;
widthwall=6*dFeO;
walldistance=100;

% graphene
c=3.2;

% water
ee=3.15;

% box size
txlo=-20*dFeO; txhi=-txlo;
tylo=-8.5223; tyhi=-tylo;
tzlo=-walldistance/2-widthwall-5; tzhi=-tzlo;

% wall 1 FeO
PnmpA=load('../../../shared/FeO/Positions.dat');
x=txlo+dFeO/2;
y=tylo+dFeO/2;
z=-widthwall/2;
cptmol=cptmol+1;
while z<=widthwall/2-2*dFeO
	while y < tyhi
		while x < txhi
			for ii=1:length(PnmpA(:,1))
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
			end
			x=x+2*dFeO;
		end
		y=y+2*dFeO;
		x=txlo+dFeO/2;
	end
	z=z+2*dFeO;
	y=tylo+dFeO/2;
end
A(:,7)=A(:,7)+walldistance/2+widthwall/2;

% wall 2 FeO
cptmol=cptmol+1;
At=A;
At(:,7)=-At(:,7);
At(:,2)=cptmol;
At(:,1)=At(:,1)+cptatom;
A=[A; At];
cptatom=max(A(:,1));

% load  graphene
PnmpA=load('../../../shared/GrapheneParticle/Positions.dat');
PnmpA(:,5)=PnmpA(:,5)-mean(PnmpA(:,5));
PnmpA(:,6)=PnmpA(:,6)-mean(PnmpA(:,6));
PnmpA(:,7)=PnmpA(:,7)-mean(PnmpA(:,7));
PnmpA(:,3) = PnmpA(:,3)+4; % shift atom ids
BnmpA=load('../../../shared/GrapheneParticle/Bonds.dat');
AnmpA=load('../../../shared/GrapheneParticle/Angles.dat');
DnmpA=load('../../../shared/GrapheneParticle/Dihedrals.dat');
InmpA=load('../../../shared/GrapheneParticle/Impropers.dat');
% place graphene layers
x=0;
y=0;
for z=-c:c:z
	cptatom0=cptatom;
	cptmol=cptmol+1;
	for ii=1:length(BnmpA(:,1))
		cptbond=cptbond+1;
		B(cptbond,:)=[cptbond BnmpA(ii,2) BnmpA(ii,3)+cptatom BnmpA(ii,4)+cptatom];
	end
	for ii=1:length(AnmpA(:,1))
		cptangle=cptangle+1;
		Ag(cptangle,:)=[cptangle AnmpA(ii,2) AnmpA(ii,3)+cptatom AnmpA(ii,4)+cptatom AnmpA(ii,5)+cptatom];
	end
	for ii=1:length(DnmpA(:,1))
		cptdihedrals=cptdihedrals+1;
		D(cptdihedrals,:)=[cptdihedrals DnmpA(ii,2) DnmpA(ii,3)+cptatom DnmpA(ii,4)+cptatom DnmpA(ii,5)+cptatom DnmpA(ii,6)+cptatom];
	end
	for ii=1:length(InmpA(:,1))
		cptimpropers=cptimpropers+1;
		Ig(cptimpropers,:)=[cptimpropers InmpA(ii,2) InmpA(ii,3)+cptatom InmpA(ii,4)+cptatom InmpA(ii,5)+cptatom InmpA(ii,6)+cptatom];
	end

	for ii=1:length(PnmpA(:,1))
		cptatom=cptatom+1;
		A(cptatom,:)=[PnmpA(ii,1)+cptatom0 cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
	end
end

% load water
lengthAbeforewater=length(A(:,1));
clear PnmpA;
clear BnmpA;
clear AnmpA;
PnmpA=load('../../../shared/H2O_SPCE/position.dat');
PnmpA(:,3) = PnmpA(:,3)+2; % shift atom ids
BnmpA=load('../../../shared/H2O_SPCE/bond.dat');
BnmpA(:,2) = BnmpA(:,2)+2; % shift bond ids
AnmpA=load('../../../shared/H2O_SPCE/angle.dat');
AnmpA(:,2) = AnmpA(:,2)+2; % shift angle ids
% place water molecules
dx=ee; dy=ee; dz=ee;
x=txlo+dx/2;
y=tylo+dy/2;
z=-walldistance/2+ee;
while z<walldistance/2-ee
	while y < tyhi-dy/2
		while x < txhi-dz/2

			xeff=x+(rand-0.5)/10;
			yeff=y+(rand-0.5)/10;
			zeff=z+(rand-0.5)/10;
		
			for ii=1:length(PnmpA(:,1))
				Glyc(ii,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+xeff PnmpA(ii,6)+yeff PnmpA(ii,7)+zeff];
			end
			mind=1e5;
			for ii=1:length(Glyc(:,1))
				for jj=1:lengthAbeforewater
					d=sqrt((Glyc(ii,5)-A(jj,5))^2+(Glyc(ii,6)-A(jj,6))^2+(Glyc(ii,7)-A(jj,7))^2);
					if d<mind
						mind=d;
					end
				end
			end
			if mind>3.4
				cptmol=cptmol+1;
				for ii=1:length(BnmpA(:,1))
					cptbond=cptbond+1;
					B(cptbond,:)=[cptbond BnmpA(ii,2) BnmpA(ii,3)+cptatom BnmpA(ii,4)+cptatom];
				end

				for ii=1:length(AnmpA(:,1))
					cptangle=cptangle+1;
					Ag(cptangle,:)=[cptangle AnmpA(ii,2) AnmpA(ii,3)+cptatom AnmpA(ii,4)+cptatom AnmpA(ii,5)+cptatom];
				end
				for ii=1:length(PnmpA(:,1))
					cptatom=cptatom+1;
					A(cptatom,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+xeff PnmpA(ii,6)+yeff PnmpA(ii,7)+zeff];
				end
			end
			x=x+dx;
		end
		x=txlo+dx/2;
		y=y+dy;
	end
	y=tylo+dy/2;
	z=z+dz;	
end

% write data lammps
fid = fopen('data.lammps','wt');
fprintf(fid, '# System\n\n');
fprintf(fid, num2str(cptatom));
fprintf(fid, ' atoms\n');
fprintf(fid, num2str(cptbond));
fprintf(fid, ' bonds\n');
fprintf(fid, num2str(cptangle));
fprintf(fid, ' angles\n');
fprintf(fid, num2str(cptdihedrals));
fprintf(fid, ' dihedrals\n');
fprintf(fid, num2str(cptimpropers));
fprintf(fid, ' impropers\n\n');
fprintf(fid, num2str(Natomtypes));
fprintf(fid, ' atom types\n');
fprintf(fid, num2str(Nbondtypes));
fprintf(fid, ' bond types\n');
fprintf(fid, num2str(Nangletypes));
fprintf(fid, ' angle types\n');
fprintf(fid, num2str(Ndihedraltypes));
fprintf(fid, ' dihedral types\n');
fprintf(fid, num2str(Nimpropertypes));
fprintf(fid, ' improper types\n');
fprintf(fid, 'extra bond per atom ');
fprintf(fid, num2str(2));
fprintf(fid, '\n');
fprintf(fid, 'extra angle per atom ');
fprintf(fid, num2str(1));
fprintf(fid, '\n');
fprintf(fid, 'extra special per atom ');
fprintf(fid, num2str(2));
fprintf(fid, '\n');
fprintf(fid, num2str([txlo txhi]));
fprintf(fid, ' xlo xhi\n');
fprintf(fid, num2str([tylo tyhi]));
fprintf(fid, ' ylo yhi\n');
fprintf(fid, num2str([tzlo tzhi]));
fprintf(fid, ' zlo zhi\n\n');
fprintf(fid, 'Atoms\n\n');
for ii=1:length(A(:,1))
	for jj=1:7
		fprintf(fid, num2str(A(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Bonds\n\n');
for ii=1:length(B(:,1))
	for jj=1:4
		fprintf(fid, num2str(B(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Angles\n\n');
for ii=1:length(Ag(:,1))
	for jj=1:5
		fprintf(fid, num2str(Ag(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Dihedrals\n\n');
for ii=1:length(D(:,1))
	for jj=1:6
		fprintf(fid, num2str(D(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Impropers\n\n');
for ii=1:length(Ig(:,1))
	for jj=1:6
		fprintf(fid, num2str(Ig(ii,jj)));
		fprintf(fid, ' ');
	end
	fprintf(fid, '\n');
end
disp('Done writing the data file')
fclose(fid);
