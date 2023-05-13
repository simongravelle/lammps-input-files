% this file generate the PARM.lammps file

% Water
MASSwater=load('H2O_TIP4P2005/file.mass'); 
PAIRCOEFFwater=load('H2O_TIP4P2005/file.paircoeff');
BONDwater=load('H2O_TIP4P2005/file.bond');
ANGLEwater=load('H2O_TIP4P2005/file.angle');

% Toluene
MASSToluene=load('toluene/file.mass');
MASSToluene(:,1)=MASSToluene(:,1)+length(MASSwater(:,1));
PAIRCOEFFToluene=load('toluene/file.paircoeff');
PAIRCOEFFToluene(:,[1:2])=PAIRCOEFFToluene(:,[1:2])+length(MASSwater(:,1));
BONDToluene=load('toluene/file.bond');
BONDToluene(:,1)=BONDToluene(:,1)+length(BONDwater(:,1));
ANGLEToluene=load('toluene/file.angle');
ANGLEToluene(:,1)=ANGLEToluene(:,1)+length(ANGLEwater(:,1));
DIHEDRALToluene=load('toluene/file.dihedral');
IMPROPERToluene=load('toluene/file.improper');

MASS=[MASSwater; MASSToluene];
PAIRCOEFF=[PAIRCOEFFwater; PAIRCOEFFToluene];
BOND=[BONDwater; BONDToluene];
ANGLE=[ANGLEwater; ANGLEToluene];

Natomtypes=length(MASS(:,1)); 
Nbondtypes=length(BOND(:,1)); 
Nangletypes=length(ANGLE(:,1)); 
Ndihedraltypes=length(DIHEDRALToluene(:,1));
Nimpropertypes=length(IMPROPERToluene(:,1)); 

disp('Start writing the data file')
fid = fopen('PARM.lammps','wt');
fprintf(fid, '#Mass\n\n');
for ii=1:length(MASS(:,1))
		fprintf(fid, 'mass ');
		fprintf(fid, num2str(MASS(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(MASS(ii,2)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Pair Coeff\n\n');
for ii=1:length(PAIRCOEFF(:,1))
		fprintf(fid, 'pair_coeff ');
		fprintf(fid, num2str(PAIRCOEFF(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,2)));
		%fprintf(fid, ' lj/cut/coul/long ');
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,3)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,4)));

	fprintf(fid, '\n');
end

fprintf(fid, '\n');
fprintf(fid, '#Bond\n\n');
for ii=1:length(BOND(:,1))
		fprintf(fid, 'bond_coeff ');
		fprintf(fid, num2str(BOND(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(BOND(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(BOND(ii,3)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Angle\n\n');
for ii=1:length(ANGLE(:,1))
		fprintf(fid, 'angle_coeff ');
		fprintf(fid, num2str(ANGLE(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(ANGLE(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(ANGLE(ii,3)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Dihedral\n\n');
for ii=1:length(DIHEDRALToluene(:,1))
		fprintf(fid, 'dihedral_coeff ');
		fprintf(fid, num2str(DIHEDRALToluene(ii,1)));
		fprintf(fid, ' harmonic ');
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRALToluene(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRALToluene(ii,3)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRALToluene(ii,4)));
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Improper\n\n');
for ii=1:length(IMPROPERToluene(:,1))
		fprintf(fid, 'improper_coeff ');
		fprintf(fid, num2str(IMPROPERToluene(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(IMPROPERToluene(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(IMPROPERToluene(ii,3)));

	fprintf(fid, '\n');
end
disp('Done writing the data file')
fclose(fid);




