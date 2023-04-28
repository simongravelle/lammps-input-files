clear all
close all
beep off

% read the datafile
fid = fopen('../../S0nmp/data.npt.lammps');
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline([end-5:end])=[];
Natoms=str2num(tline);

tline = fgetl(fid);
tline([end-10:end])=[];
Tatoms=str2num(tline);

tline = fgetl(fid);
tline([end-5:end])=[];
Nbonds=str2num(tline);

tline = fgetl(fid);
tline([end-10:end])=[];
Tbonds=str2num(tline);

tline = fgetl(fid);
tline([end-6:end])=[];
Nangles=str2num(tline);

tline = fgetl(fid);
tline([end-11:end])=[];	
Tangles=str2num(tline);

tline = fgetl(fid);
tline([end-9:end])=[];
Ndihedrals=str2num(tline);

tline = fgetl(fid);
tline([end-14:end])=[];
Tdihedrals=str2num(tline);

tline = fgetl(fid);
tline([end-9:end])=[];
Nimpropers =str2num(tline);

tline = fgetl(fid);
tline([end-14:end])=[];
Timpropers=str2num(tline);

tline = fgetl(fid);

tline = fgetl(fid);
tline([end-7:end])=[];
xcoor=str2num(tline); Lx=xcoor(2)-xcoor(1);
tline = fgetl(fid);
tline([end-7:end])=[];
ycoor=str2num(tline); Ly=ycoor(2)-ycoor(1);
tline = fgetl(fid);
tline([end-7:end])=[];
zcoor=str2num(tline); Lz=zcoor(2)-zcoor(1);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
%
%for ii=1:Tatoms
%	tline = fgetl(fid);
%	%tline([end-4:end])=[];
%	MassCoeff(ii,:)=str2num(tline);
%end
%
%tline = fgetl(fid);
%tline = fgetl(fid);
%tline = fgetl(fid);
%
for ii=1:Natoms
	tline = fgetl(fid);
	%tline([end-8:end])=[];
	Positions(ii,:)=str2num(tline);
end
%
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
for ii=1:Natoms
	tline = fgetl(fid);
	%tline([end-8:end])=[];
	Vel(ii,:)=str2num(tline);
end
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
%
for ii=1:Nbonds
	tline = fgetl(fid);
	Bonds(ii,:)=str2num(tline);
end
%
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
%
for ii=1:Nangles
	tline = fgetl(fid);
	Angles(ii,:)=str2num(tline);
end
%
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
%
for ii=1:Ndihedrals
	tline = fgetl(fid);
	Dihedrals(ii,:)=str2num(tline);
end
%
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
%
for ii=1:Nimpropers
	tline = fgetl(fid);
	Impropers(ii,:)=str2num(tline);
end


coor=[xcoor; ycoor; zcoor];

fid2 = fopen('Coordinate.dat','wt');
for ii=1:length(coor(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(coor(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(coor(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);

fid2 = fopen('Positions.dat','wt');
for ii=1:length(Positions(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==3
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==4
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==5
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==6
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==7
			fprintf(fid2, num2str(Positions(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);

fid2 = fopen('Bonds.dat','wt');
for ii=1:length(Bonds(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(ii));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(Bonds(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==3
			fprintf(fid2, num2str(Bonds(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==4
			fprintf(fid2, num2str(Bonds(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);

fid2 = fopen('Angles.dat','wt');
for ii=1:length(Angles(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(ii));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(Angles(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==3
			fprintf(fid2, num2str(Angles(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==4
			fprintf(fid2, num2str(Angles(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==5
			fprintf(fid2, num2str(Angles(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);

fid2 = fopen('Dihedrals.dat','wt');
for ii=1:length(Dihedrals(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(ii));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(Dihedrals(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==3
			fprintf(fid2, num2str(Dihedrals(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==4
			fprintf(fid2, num2str(Dihedrals(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==5
			fprintf(fid2, num2str(Dihedrals(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==6
			fprintf(fid2, num2str(Dihedrals(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);

fid2 = fopen('Impropers.dat','wt');
for ii=1:length(Impropers(:,1))
	for jj=1:7
		if jj==1
			fprintf(fid2, num2str(ii));
			fprintf(fid2, '	');
		elseif jj==2
			fprintf(fid2, num2str(Impropers(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==3
			fprintf(fid2, num2str(Impropers(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==4
			fprintf(fid2, num2str(Impropers(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==5
			fprintf(fid2, num2str(Impropers(ii,jj)));
			fprintf(fid2, '	');
		elseif jj==6
			fprintf(fid2, num2str(Impropers(ii,jj)));
			fprintf(fid2, '	');
		end
	end
	fprintf(fid2, '\n');
end
fclose(fid2);


