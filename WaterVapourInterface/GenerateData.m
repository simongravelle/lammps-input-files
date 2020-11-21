% this octave script generate the file data.lammps to be read from lammps input script. it contains a slack of water molecule

clear all
close all

cptmol=0;
cptatom=0;
cptbond=0;
cptangle=0;

txlo=-10; txhi=-txlo;
tylo=-10; tyhi=-tylo;
tzlo=-30; tzhi=-tzlo;

dx=3.4;

PnmpA=load('./WaterMolecule/Position.dat');
BnmpA=load('./WaterMolecule/Bond.dat');
AnmpA=load('./WaterMolecule/Angle.dat');

x=txlo+dx/4;
y=tylo+dx/4;
z=tzlo/2;
while z<=tzhi/2
	while y < tyhi-dx/2
		while x < txhi-dx/2
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
			x=x+dx;
		end
		x=txlo+dx/4;
		y=y+dx;
	end
	y=tylo+dx/4;
	z=z+dx;	
end


fid = fopen('data.lammps','wt');
fprintf(fid, '# System\n\n');
fprintf(fid, num2str(cptatom));
fprintf(fid, ' atoms\n');
fprintf(fid, num2str(cptbond));
fprintf(fid, ' bonds\n');
fprintf(fid, num2str(cptangle));
fprintf(fid, ' angles\n');
fprintf(fid, num2str(2));
fprintf(fid, ' atom types\n');
fprintf(fid, num2str(1));
fprintf(fid, ' bond types\n');
fprintf(fid, num2str(1));
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
fclose(fid);



