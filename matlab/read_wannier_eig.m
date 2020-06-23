function E = read_wannier_eig(filename,Nb,Nk)
% function to read data from wannier90 eig file
% CALL: E = READ_WANNIER_EIG('wannier90.eig',Nb,Nk)
%
% input args:
% - filename: name of the wannier90 eig file, default: 'wannier90.eig'
% - Nb: number of bands
% - Nk: number of k-points
%
% output args:
% - E: energy matrix, size [Nb Nk]

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


fid = fopen(filename);
E = zeros(Nb,Nk);

for ck=1:Nk
    for cb=1:Nb
        dat = sscanf(fgets(fid),'%i %i %lf');
        if cb~=dat(1) || ck~=dat(2)
            error('gopfertami')
        end
        E(cb,ck) = dat(3);
    end
end

end
