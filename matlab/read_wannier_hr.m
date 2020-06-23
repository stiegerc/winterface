function [H,R,Nw,NR] = read_wannier_hr(filename)
% function to read hamiltonian data from wannier90 output
% CALL: [H,R,Nw,NR] = READ_WANNIER_HR('wannier90_hr.dat')
%
% input args:
% - filename: name of the wannier90 hr file, default: 'wannier90_hr.dat'
%
% output args:
% - H: hamiltonian cellarray in blocks H{R}
% - R: R vectors for each block in H
% - Nw: number of wannier functions
% - NR: number of R vectors

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


% default argument
if nargin<1
    filename = 'wannier90_hr.dat';
end


fid = fopen(filename);

fgets(fid);
Nw = sscanf(fgets(fid),'%lf');
NR = sscanf(fgets(fid),'%lf');

for skip=1:ceil(NR/15)
    fgets(fid);
end

R = zeros(3,NR);
H = cell(1,NR);

buff = fgets(fid);
for cr=1:NR
    R(:,cr) = sscanf(buff,'%lf %lf %lf',3);
    
    H{cr} = zeros(Nw) + 1i*zeros(Nw);
    for cwn=1:Nw
        for cwm=1:Nw
            hdat = sscanf(buff,'%*i %*i %*i %*i %*i %lf %lf',2);
            H{cr}(cwm,cwn) = hdat(1) + 1i*hdat(2);
            buff = fgets(fid);
        end
    end
end

fclose(fid);
