function [B,Ap,id,nn,bl] = read_lattice_dat(filename)
% function to read data from omen style lattice_dat file
% CALL: [B,Ap,id,nn,bl] = READ_LATTICE_DAT(filename)
%
% input args:
% - filename: name of the lattice_dat wout file, default: 'lattice_dat'
%
% output args:
% - B: lattice vectors as columns, size 3x3
% - Ap: atomic positions in cartesian as columns, size 3xNa
% - id: atomic type strings, size 1xNa
% - nn: number of next neighbors
% - bl: bond length

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'lattice_dat';
end


%% open file
fid = fopen(filename);


%% read Na, nn
ibuff = sscanf(fgets(fid),'%lf %lf',2);
Na = ibuff(1); nn=ibuff(2);
fgets(fid);


%% read B
B = zeros(3);
B(:,1) = sscanf(fgets(fid),'%lf %lf %lf',3);
B(:,2) = sscanf(fgets(fid),'%lf %lf %lf',3);
B(:,3) = sscanf(fgets(fid),'%lf %lf %lf',3);
fgets(fid);


%% read Ap and id
Ap = zeros(3,Na);
id = cell(1,Na);
for ca=1:Na
    buff = fgets(fid);
    id{ca} = sscanf(buff,'%s',1);
    Ap(:,ca) = sscanf(buff,'%*s %lf %lf %lf',3);
end
fgets(fid);


%% read bl
bl = sscanf(fgets(fid),'%lf',1);


%% close file
fclose(fid);

end
