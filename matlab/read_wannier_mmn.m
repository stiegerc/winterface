function [MMN,mm,Nb,Nk,nntot] = read_wannier_mmn(filename)
% function to read data from wannier90 mmn file
% CALL: [MMN,mm,Nb,Nk,nntot] = READ_WANNIER_MMN(filename)
%
% input args:
% - filename: name of the wannier90 mmn file, default: 'wannier90.mmn'
%
% output args:
% - MMN: cell array containing overlaps, size [1 Nk*nntot]
%        each element size [Nb Nb]
% - mm: meta matrix, containing [ik1 ik2 G] on each row, size [Nk*nntot 5]
% - Nb: # bands
% - Nk: # kpoints
% - nntot: # next neighbors

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'wannier90.mmn';
end

%% open file
fid = fopen(filename);

%% read Nb,Nk,nntot
fgets(fid);
dat = sscanf(fgets(fid),'%i %i %i',3);
Nb = dat(1);
Nk = dat(2);
nntot = dat(3);

%% allocate
MMN = cell(1,Nk*nntot);
mm = zeros(Nk*nntot,5);

%% read data
for c=1:Nk*nntot
    MMN{c} = zeros(Nb,Nb);
    mm(c,:) = sscanf(fgets(fid),'%i %i %i %i %i',5);
    for cb=1:Nb*Nb
        dat = sscanf(fgets(fid),'%lf %lf',2);
        MMN{c}(cb) = dat(1)+1i*dat(2);
    end
end

%% close file
fclose(fid);

end
