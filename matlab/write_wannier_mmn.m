function write_wannier_mmn(MMN,mm,Nb,Nk,nntot,filename)
% function to write wannier90 mmn data
% CALL: WRITE_WANNIER_MMN(MMN,mm,Nb,Nk,nntot,filename)
%
% input args:
% - MMN: MMN data, cell of size [1 Nk*nntot], each cell [Nb Nb]
% - mm: meta matrix, containing [ik1 ik2 G] on each row, size [Nk*nntot 5]
% - Nb: # bands
% - Nk: # kpoints
% - nntot: # next neighbors
% - filename: name of the wannier90 amn file, default: 'wannier90.mmn'

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<6
    filename = 'wannier90.mmn';
end

%% bitching
if numel(MMN)~=Nk*nntot
    error('MMN needs to be of size [1 Nk*nntot]')
end
for c=1:numel(MMN)
    if size(MMN{c})~=[Nb Nb]
        error('MMN contains bad data block of bad size')
    end
end
if size(mm)~=[Nk*nntot 5]
    error('mm bad size')
end

%% write to file
fid = fopen(filename,'w');

fprintf(fid,[datestr(datetime('now')) '\n']);
fprintf(fid,'%12.0f %12.0f %12.0f\n',Nb,Nk,nntot);

for c=1:Nk*nntot
    fprintf(fid,'%5.0f %5.0f %5.0f %5.0f %5.0f\n', ...
        mm(c,1),mm(c,2),mm(c,3),mm(c,4),mm(c,5));
    for d=1:Nb*Nb
        fprintf(fid,'%18.12f %18.12f\n', ...
            real(MMN{c}(d)),imag(MMN{c}(d)));
    end
end

end
