function [AMN,Nk,Nb,Norb] = read_wannier_amn(filename)
% function to read data from wannier90 amn file
% CALL: [AMN,Nk,Nb,Norb] = READ_WANNIER_AMN('wannier90.amn')
%
% input args:
% - filename: name of the wannier90 amn file, default: 'wannier90.amn'
%
% output args:
% - AMN: array of size [Nk Nb Norb]
% - Nk: # kpoints
% - Nb: # bands
% - Norb: # orbitals

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'wannier90.amn';
end

%% open file
fid = fopen(filename);

%% read Nb,Nk,nntot
fgets(fid);
dat = sscanf(fgets(fid),'%i %i %i',3);
Nb = dat(1);
Nk = dat(2);
Norb = dat(3);

%% allocate and read data
AMN = zeros(Nk,Nb,Norb)+1i*zeros(Nk,Nb,Norb);
for ck=1:Nk
    for corb=1:Norb
        for cb=1:Nb
            dat = sscanf(fgets(fid),'%*i %*i %*i %lf %lf',2);
            AMN(ck,cb,corb) = dat(1)+1i*dat(2);
        end
    end
end
