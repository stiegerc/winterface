function [Ef,E,k,kpl,occ] = read_outcar(filename,kI,bI)
% function to read the fermi energy, E(k) and occupations from OUTCAR
% CALL: [Ef,E,k,kpl,occ] = READ_OUTCAR('OUTCAR')
%
% input args:
% - filename: name of the OUTCAR file, default: 'OUTCAR'
% - kI: k point indices to be included, default: 1:Nk
% - bI: band indices to be included, default: 1:Nb
%
% output args:
% - Ef: the fermi energy
% - E: energy for each band and k-point, size NbIxNkI
% - k: k-points in reciprocal basis, size 3xNkI
% - kpl: k plotting vector for easy plotting, size NbIxNkI
% - occ: occupation for each band and k-point, size NbIxNkI

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'OUTCAR';
end
if nargin<2
    kI = [];
end
if nargin<3
    bI = [];
end

%% open file
fid = fopen(filename);


%% read lattice vectors
buff=[];
while isempty(strfind(buff,'Lattice vectors'))
    buff=fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
B = zeros(3);
fgets(fid);
B(:,1) = sscanf(fgets(fid),'%*s %*c %*c %lf %*c %lf %*c %lf %*c',3);
B(:,2) = sscanf(fgets(fid),'%*s %*c %*c %lf %*c %lf %*c %lf %*c',3);
B(:,3) = sscanf(fgets(fid),'%*s %*c %*c %lf %*c %lf %*c %lf %*c',3);


%% search for NKPTS, read #k-points, #bands
buff=[];
while isempty(strfind(buff,'NKPTS'))
    buff=fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
Nk = sscanf(buff,'%*s %*s %*s %lf',1);
Nb = sscanf(buff, ...
    '%*s %*s %*s %*f %*s %*s %*s %*s %*s %*f %*s %*s %*s %*s %lf',1);


%% search for E-fermi, read fermi energy
buff=[];
while isempty(strfind(buff,'E-fermi'))
    buff=fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
Ef = sscanf(buff,'%*s %*s %lf');

if nargout==1
    return
end


%% read k-points and energies
k = zeros(3,Nk);
E = zeros(Nb,Nk);
occ = zeros(Nb,Nk);

fgets(fid);
for ck=1:Nk
    buff=[];
    while isempty(strfind(buff,'k-point'))
        buff=fgets(fid);
        if feof(fid)
            error('unexpected end of file')
        end
    end
    k(:,ck) =  sscanf(buff,'%*s %*f %*s %lf %lf %lf',3);
    
    fgets(fid);
    for cb=1:Nb
        buff = fgets(fid);
        fbuff = sscanf(buff,'%*f %lf %lf',2);
        
        E(cb,ck) = fbuff(1);
        occ(cb,ck) = fbuff(2);
    end
end

fclose(fid);


%% cut out relevant kpt and bands
if isempty(kI)
    kI = 1:Nk;
else
    if min(kI)<1 || max(kI)>Nk
        error('invalid range in kI')
    end
    if any(round(kI)~=kI)
        error('non integer found in kI')
    end
    kI = sort(kI);
end
if isempty(bI)
    bI = 1:Nb;
else
    if min(bI)<1 || max(bI)>Nb
        error('invalid range in bI')
    end
    if any(round(bI)~=bI)
        error('non integer found in bI')
    end
    bI = sort(bI);
end

k = k(:,kI);
E = E(bI,kI);
occ = occ(bI,kI);


%% determine kpl
R = 2*pi*inv(B)';
kpl = zeros(1,size(k,2));
for ck=2:size(k,2)
    kpl(ck) = kpl(ck-1)+norm(R*(k(:,ck)-k(:,ck-1)));
end
kpl = ones(size(E,1),1)*kpl;
