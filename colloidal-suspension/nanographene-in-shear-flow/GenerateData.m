clear all
close all

Natomtypes=5;
Nbondtypes=9;
Nangletypes=3;
Ndihedraltypes=2;  
Nimpropertypes=1;

cptatom=0;
cptmol=0;
cptbond=0;
cptangle=0;
cptdihedrals=0;
cptimpropers=0;

wallZ0=45;
sigma=3.374; % A
unit=2*sigma/sqrt(2);
p=sqrt(2)*unit/2*0.7071067811865476;
txlo=-unit*3; txhi=-txlo;
tylo=-unit*3; tyhi=-tylo;
tzlo=-wallZ0-15; tzhi=-tzlo;

% first walls
x=txlo+unit/4;
y=tylo+unit/4;
z=wallZ0;
cptmol=cptmol+1;

while y < tyhi
	while x < txhi
		while z <= wallZ0+unit
			x0=x;
			y0=y;
			z0=z+unit/2;
			if x0<txhi && x0>=txlo && y0<tyhi && y0>=tylo
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol 3 0 x0 y0 z0];
			end
			x0=x;
			y0=y+unit/2;
			z0=z;
			if x0<txhi && x0>=txlo && y0<tyhi && y0>=tylo
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol 3 0 x0 y0 z0];
			end
			x0=x+unit/2;
			y0=y;
			z0=z;
			if x0<txhi && x0>=txlo && y0<tyhi && y0>=tylo
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol 3 0 x0 y0 z0];
			end
			x0=x+unit/2;
			y0=y+p;
			z0=z+p;
			if x0<txhi && x0>=txlo && y0<tyhi && y0>=tylo
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol 3 0 x0 y0 z0];
			end
			z=z+unit;
		end
		x=x+unit;
		z=wallZ0;
	end
	x=txlo+unit/4;
	y=y+unit;
end

% second walls
Sym=A; Sym(:,2)=Sym(:,2)+1; cptmol=cptmol+1;
Sym(:,1)=Sym(:,1)+max(Sym(:,1));
Sym(:,7)=-Sym(:,7);
A=[A;Sym];
cptatom=max(A(:,1));

% Hexaben molecule
PnmpA=load('../../shared/Hexabenzocoronene/Position.dat');
BnmpA=load('../../shared/Hexabenzocoronene/Bond.dat');
AnmpA=load('../../shared/Hexabenzocoronene/Angle.dat');
DnmpA=load('../../shared/Hexabenzocoronene/Dihedral.dat');
InmpA=load('../../shared/Hexabenzocoronene/Impropers.dat');

x=0; y=0; z=0;
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
	A(cptatom,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
end

% water molecules
nbCarWall=length(A(:,1));
PnmpA=load('../../shared/H2O_TIP4P2005/position.dat');
PnmpA(:,3) = PnmpA(:,3) + 3;
BnmpA=load('../../shared/H2O_TIP4P2005/bond.dat');
BnmpA(:,2) = BnmpA(:,2) + 8;
AnmpA=load('../../shared/H2O_TIP4P2005/angle.dat');
AnmpA(:,2) = AnmpA(:,2) + 2;

dx=3.4;
dy=3.4;
dz=3.4;
x=txlo+dx/2;
y=tylo+dy/2;
z=-wallZ0+dz;
while z<wallZ0-dz
	while y < tyhi-dy/2
		while x < txhi-dx/2
			dmin=10;
			for cpt=1:nbCarWall
				d=sqrt((A(cpt,5)-x)^2+(A(cpt,6)-y)^2+(A(cpt,7)-z)^2);
				dmin=min(d,dmin);	
			end
			if dmin>3
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
					A(cptatom,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
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

% write data.lammps
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
fprintf(fid, '\n\n');

fprintf(fid, '# Membrane\n\n');

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
for ii=1:length(Ig(:,1))%%
	for jj=1:6%
		fprintf(fid, num2str(Ig(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end



disp('Done writing the data file')

fclose(fid);










































