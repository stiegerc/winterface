function write_wannier_amn(AMN,filename)
% function to write wannier90 amn data
% CALL: WRITE_WANNIER_AMN(AMN,filename)
%
% input args:
% - AMN: AMN data, size [Nk Nb Norb]
% - filename: name of the wannier90 amn file, default: 'wannier90.amn'

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<2
    filename = 'wannier90.amn';
end

%% bitching
if length(size(AMN))<3
    error('AMN needs to be of size [Nk Nb Norb]')
end
Nk = size(AMN,1); Nb = size(AMN,2); Norb = size(AMN,3);
if (Nb<Norb)
    disp('warning: Nb<Norb')
end

%% write to file
fid = fopen(filename,'w');

fprintf(fid,[datestr(datetime('now')) '\n']);
fprintf(fid,'%12.0f %12.0f %12.0f\n',Nb,Nk,Norb);

for ck=1:Nk
    for co=1:Norb
        for cb=1:Nb
            fprintf(fid,'%5.0f %5.0f %5.0f %18.12f %18.12f\n',...
                cb,co,ck,real(AMN(ck,cb,co)),imag(AMN(ck,cb,co)));
        end
    end
end
