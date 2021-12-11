% this Octave/Matlab algorithm add Toluene and Water molecule to a cubid box

clear all
close all

GenerateParamFile

cptmol=0;
cptatom=0;
cptbond=0;
cptangle=0;
cptdihedrals=0;
cptimpropers=0;

txlo=-30; txhi=-txlo; Lx=txhi-txlo;
tylo=-80; tyhi=-tylo; Ly=tyhi-tylo;
tzlo=-30; tzhi=-tzlo; Lz=tzhi-tzlo;

	
%%%%%%%%%%%%%%%%%%%%%
% toluene molecules %
%%%%%%%%%%%%%%%%%%%%%

PnmpA=load('../../shared/toluene/Positions.dat');
PnmpA(:,3)=PnmpA(:,3)+2; % shift to account for water molecule (O and H) 
PnmpA(:,5)=PnmpA(:,5)-mean(PnmpA(:,5)); % recenter the Toluene molecule
PnmpA(:,6)=PnmpA(:,6)-mean(PnmpA(:,6)); % recenter the Toluene molecule
PnmpA(:,7)=PnmpA(:,7)-mean(PnmpA(:,7)); % recenter the Toluene molecule
BnmpA=load('../../shared/toluene/Bonds.dat');
BnmpA(:,2)=BnmpA(:,2)+1; % shift OH water
AnmpA=load('../../shared/toluene/Angles.dat');
AnmpA(:,2)=AnmpA(:,2)+1; % shift HOH water
DnmpA=load('../../shared/toluene/Dihedrals.dat');
DnmpA(:,2)=DnmpA(:,2);
InmpA=load('../../shared/toluene/Impropers.dat');
InmpA(:,2)=InmpA(:,2);

dx=8; dy=6; dz=7;

% add Toluene molecule at random position if no intersection with other molecule
Ntrial=0;
failinsert=0;
while Ntrial<1e4
	x=rand*(txhi-txlo)+txlo;
	y=(rand*(tyhi-tylo)+tylo)*0.9;
	z=rand*(tzhi-tzlo)+tzlo;
	Ntrial=Ntrial+1;
	
	mind=1e5;
	TolMol=PnmpA;
	TolMol(:,5)=TolMol(:,5)+x;
	TolMol(:,6)=TolMol(:,6)+y;
	TolMol(:,7)=TolMol(:,7)+z;
	if cptmol>0
		for ii=1:length(TolMol(:,1))
			for jj=1:length(A(:,1))
				for x0=-Lx:Lx:Lx
					for y0=-Ly:Ly:Ly
						for z0=-Lz:Lz:Lz
							d=sqrt((x0+TolMol(ii,5)-A(jj,5))^2+(y0+TolMol(ii,6)-A(jj,6))^2+(z0+TolMol(ii,7)-A(jj,7))^2);
							if d<mind
								mind=d;
							end
						end
					end
				end
			end
		end
	end

	if mind>3.4 % no interection -> add the water molecule
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
	else
		failinsert=failinsert+1;
	end
	
	% print ou the success/fail of water insertion
	if mod(Ntrial,100)==0
		X=['toluene trial n°', num2str(Ntrial),', success = ', num2str(cptmol),', fail = ', num2str(failinsert)];
		disp(X);
	end
end

X=[num2str(cptmol),' toluene molecule have been added'];
disp(X)


%%%%%%%%%%%%%%%%%%%
% water molecules %
%%%%%%%%%%%%%%%%%%%

clear PnmpA; 
clear BnmpA;
clear AnmpA;
clear DnmpA;
clear InmpA;

PnmpA=load('../../shared/H2O_TIP4P2005/Position.dat');
BnmpA=load('../../shared/H2O_TIP4P2005/Bond.dat');
AnmpA=load('../../shared/H2O_TIP4P2005/Angle.dat');

% add water molecule at random position if no intersection with other molecule
Ntrial=0;
cptwat=0;
failinsert=0;
while Ntrial<5e4
	%rho=rand*(R+3.4);
	%theta=rand*2*pi;
	x=rand*(txhi-txlo)+txlo;
	y=(rand*(tyhi-tylo)+tylo)*0.9;
	z=rand*(tzhi-tzlo)+tzlo;
	Ntrial=Ntrial+1;
	
	mind=1e5;
	for jj=1:length(A(:,1))
		for x0=-Lx:Lx:Lx
			for y0=-Ly:Ly:Ly
				for z0=-Lz:Lz:Lz
					d=sqrt((x0+x-A(jj,5))^2+(y0+y-A(jj,6))^2+(z0+z-A(jj,7))^2);
					if d<mind
						mind=d;
					end
				end
			end
		end
	end

	if mind>3.4 % no interection -> add the water molecule
		cptmol=cptmol+1;
		cptwat=cptwat+1;
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
	else
		failinsert=failinsert+1;
	end
	
	if mod(Ntrial,100)==0
		% print ou the success/fail of water insertion
		X=['water trial n°', num2str(Ntrial),', success = ', num2str(cptwat),', fail = ', num2str(failinsert)];
		disp(X);
	end
end

X=[num2str(cptwat),' water molecule have been added'];
disp(X)

% generate the lammps data file
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
fclose(fid);





























