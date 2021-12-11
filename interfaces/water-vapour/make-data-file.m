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

Ph2O=load('../../shared/H2O_TIP4P2005/position.dat');
Bh2O=load('../../shared/H2O_TIP4P2005/bond.dat');
Ah2O=load('../../shared/H2O_TIP4P2005/angle.dat');

x=txlo+dx/4;
y=tylo+dx/4;
z=tzlo/2;
while z<=tzhi/2
	while y < tyhi-dx/2
		while x < txhi-dx/2
			cptmol=cptmol+1; 
			for ii=1:length(Bh2O(:,1))
				cptbond=cptbond+1;
				B(cptbond,:)=[cptbond Bh2O(ii,2) Bh2O(ii,3)+cptatom Bh2O(ii,4)+cptatom];
			end
			for ii=1:length(Ah2O(:,1))
				cptangle=cptangle+1;
				Ag(cptangle,:)=[cptangle Ah2O(ii,2) Ah2O(ii,3)+cptatom Ah2O(ii,4)+cptatom Ah2O(ii,5)+cptatom];
			end
			for ii=1:length(Ph2O(:,1))
				cptatom=cptatom+1;
				A(cptatom,:)=[cptatom cptmol Ph2O(ii,3) Ph2O(ii,4) Ph2O(ii,5)+x Ph2O(ii,6)+y Ph2O(ii,7)+z];
			end
			x=x+dx;
		end
		x=txlo+dx/4;
		y=y+dx;
	end
	y=tylo+dx/4;
	z=z+dx;	
end

X = ['The number of water molecule is ',num2str(cptmol)];
disp(X)

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



