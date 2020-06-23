function [U,Udis,k] = read_wannier_umats(seedname)
% function to read wannier90 U and Udis matrices
% CALL: [U,Udis,k] = READ_WANNIER_UMATS('wannier90')
%
% input args:
% - seedname: wannier90 seed, default 'wannier90'
%
% output args:
% - U: cell array, size [1 Nk]
% - Udis: cell array, size [1 Nk] or [0] if no disentanglement
% - k: k vectors, size [3 Nk]

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default arguments
if nargin<1
    seedname = 'wannier90';
end

%% read U
fid = fopen([seedname '_u.mat']);

fgets(fid); % comment
N = sscanf(fgets(fid),'%i %i %i',3); % meta nfo

U = cell(1,N(1));
k = zeros(3,N(1));
for ck=1:N(1);
    fgets(fid);
    
    k(:,ck) = sscanf(fgets(fid),'%lf %lf %lf');
    U{ck} = zeros(N(3),N(2))+1i*zeros(N(3),N(2));
    
    for cu=1:N(3)*N(2)
        dat = sscanf(fgets(fid),'%lf %lf',2);
        U{ck}(cu) = dat(1)+1i*dat(2);
    end
end

%% read Udis
if exist([seedname '_u_dis.mat'],'file')==2
    fid = fopen([seedname '_u_dis.mat']);
    fgets(fid); %comment
    
    N = sscanf(fgets(fid),'%i %i %i',3); % meta nfo
    
    Udis = cell(1,N(1));
    for ck=1:N(1);
        fgets(fid);

        k_ = sscanf(fgets(fid),'%lf %lf %lf');
        if (any(k_-k(:,ck))>1e-6)
            error('k mismatch in Udis file')
        end
        
        Udis{ck} = zeros(N(3),N(2))+1i*zeros(N(3),N(2));

        for cu=1:N(3)*N(2)
            dat = sscanf(fgets(fid),'%lf %lf',2);
            Udis{ck}(cu) = dat(1)+1i*dat(2);
        end
    end
else
    Udis=cell(0);
end

end
