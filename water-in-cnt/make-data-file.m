% this Octave/Matlab algorithm add Water molecule in a CNT

clear all
close all

Natomtypes=3;
Nbondtypes=1;
Nangletypes=1;

% initialise counter
cptmol=0;
cptatom=0;
cptbond=0;
cptangle=0;

% set box size
txlo=-10; txhi=-txlo; Lx=txhi-txlo;
tylo=-10; tyhi=-tylo; Ly=tyhi-tylo;
tzlo=-27.8696-0.89990265356/2; tzhi=-tzlo; Lz=tzhi-tzlo;

% import CNT atoms' positions
A=load('CNT/position.dat');
A(:,3)=A(:,3)+2; % shift to account for water molecule (O and H) 
A(:,5)=A(:,5)-mean(A(:,5)); % recenter the CNT
A(:,6)=A(:,6)-mean(A(:,6)); % recenter the CNT
A(:,7)=A(:,7)-mean(A(:,7)); % recenter the CNT
cptatom=length(A(:,1));
cptmol=cptmol+1;
A(:,2)=cptmol;

% import water molecules details
Ph2o=load('H2O_TIP4P2005/position.dat');
Bh2o=load('H2O_TIP4P2005/bond.dat');
Ah2o=load('H2O_TIP4P2005/angle.dat');

% add a few water molecules inside the CNT
x=0;
y=0;
for z=tzlo+3.4/2:3.4:tzhi-3.4/2
	cptmol=cptmol+1;
	for ii=1:length(Bh2o(:,1))
		cptbond=cptbond+1;
		B(cptbond,:)=[cptbond Bh2o(ii,2) Bh2o(ii,3)+cptatom Bh2o(ii,4)+cptatom];
	end
	for ii=1:length(Ah2o(:,1))
		cptangle=cptangle+1;
		Ag(cptangle,:)=[cptangle Ah2o(ii,2) Ah2o(ii,3)+cptatom Ah2o(ii,4)+cptatom Ah2o(ii,5)+cptatom];
	end
	for ii=1:length(Ph2o(:,1))
		cptatom=cptatom+1;
		A(cptatom,:)=[cptatom cptmol Ph2o(ii,3) Ph2o(ii,4) Ph2o(ii,5)+x Ph2o(ii,6)+y Ph2o(ii,7)+z];
	end
end

X = ['The number of water molecule is ',num2str(cptmol)];
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
fprintf(fid, num2str(Natomtypes));
fprintf(fid, ' atom types\n');
fprintf(fid, num2str(Nbondtypes));
fprintf(fid, ' bond types\n');
fprintf(fid, num2str(Nangletypes));
fprintf(fid, ' angle types\n');
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
fclose(fid);
