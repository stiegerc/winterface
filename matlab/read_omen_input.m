function [B,Ap,T,mass,id,Norb,nn,bl,Eg,cb,vb] = read_omen_input(mp_file,ld_file)
% function to read data from omen style mat_par and lattice_dat files
% CALL: [B,Ap,T,mass,id,Norb,nn,bl,Eg,cb,vb] = read_omen_input(mp_file,ld_file)
%
% input args:
% - mp_file: name of the mat_par file, default: 'ph_mat_par'
% - ld_file: name of the lattice_dat file: 'lattice_dat'
%
% output args:
% - B: lattice vectors as columns, size 3x3
% - Ap: atomic positions in cartesian as columns, size 3xNa
% - T: type vector corresponding to Ap, size 1xNa
% - mass: vector of masses, size 1xnumel(unique(T))
% - id: cell array of type strings, size 1xnumel(unique(T))
% - Norb: number of orbitals for each type, size 1xNa
% - nn: number of next neighbors
% - bl: bond length
% - Eg: energy gap on eV
% - cb: conduction bands edge in eV
% - vb: valence band edge in eV

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    mp_file = 'ph_mat_par';
end
if nargin<2
    ld_file = 'lattice_dat';
end


%% read mat_par
fid = fopen(mp_file);

% number of types
NT = sum(sscanf(fgets(fid),'%lf',2));

% energies
E = sscanf(fgets(fid),'%lf',3);
Eg = E(1); cb = E(2); vb = E(3);

% Norb
Norb = [];
for c=1:ceil(NT/20)
    Norb = [Norb sscanf(fgets(fid),'%lf',20)'];
end

% mass
mass = sscanf(fgets(fid),'%lf',NT);

fclose(fid);


%% read lattice_dat
fid = fopen(ld_file);

% read Na, nn
ibuff = sscanf(fgets(fid),'%lf %lf',2);
Na = ibuff(1); nn=ibuff(2);
fgets(fid);

% read B
B = zeros(3);
B(:,1) = sscanf(fgets(fid),'%lf %lf %lf',3);
B(:,2) = sscanf(fgets(fid),'%lf %lf %lf',3);
B(:,3) = sscanf(fgets(fid),'%lf %lf %lf',3);
fgets(fid);

% read Ap and all id
Ap = zeros(3,Na); id_ = cell(1,Na);
for ca=1:Na
    buff = fgets(fid);
    id_{ca} = sscanf(buff,'%s',1);
    Ap(:,ca) = sscanf(buff,'%*s %lf %lf %lf',3);
end

% get unique id in order of occurrence
[id,ia] = unique(id_);
[~,ic] = sort(ia);
id = id(ic);

% get T using id
T = zeros(1,Na);
for ca=1:Na
    T(ca) = find(strcmp(id,id_{ca}));
end

% read bl
fgets(fid);
bl = sscanf(fgets(fid),'%lf',1);

fclose(fid);

end
