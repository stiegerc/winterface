function [k,B,w,Nk,Nb,Ni] = read_procar(filename)
% function to read PROCAR files from VASP
% CALL: [k,B,w,Nk,Nb,Ni] = READ_PROCAR(filename)
%
% input args:
% - filename: filename of the PROCAR file, default 'PROCAR'
%
% output args:
% - k: kpoints, size 3xNk
% - B: band data, array of length Nk, each element holds:
%   - E: energies in eV for each k-point, size 1xNb
%   - occ: occupation of state in E, size 1xNb
%   - s,p,d,...: orbital projection coefficients, size NixNb
% - w: weights, a weight factor for each kpoint, size 1xNk
% - Nk: # of kpoints in PROCAR
% - Nb: # of bands in PROCAR
% - Ni: # of ions in PROCAR

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'PROCAR';
end

%% open file, define format switches
fid = fopen(filename);
buff = fgets(fid);

% check for decomposed orbitals
lm_decomp = ~isempty(strfind(buff,'lm decomposed'));

% check if we have f orbitals
fgets(fid);
buff=[];
while isempty(strfind(buff,'ion'))
    buff=fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
have_f = ~isempty(strfind(buff,'f'));

% check of we have charge entries
buff=[];
while isempty(strfind(buff,'charge')) && isempty(strfind(buff,'band'))
    buff=fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
have_charge = ~isempty(strfind(buff,'charge'));

fclose(fid);

%% define format
format = '%lf %lf'; % always have s, tot
if ~lm_decomp
    format = [format ' %lf %lf']; % add p,d
    if have_f
        format = [format ' %lf']; % add f
    end
else
    format = [format ' %lf %lf %lf %lf %lf %lf %lf %lf']; % add py,pz,px,dxy,dyz,dz2,dxz,x2-y2
    if have_f
        format = [format ' %lf %lf %lf %lf %lf %lf %lf']; % add fy3x2,fxyz,fyz2,fz3,fxz2,fzx2,fx3
    end
end

%% reopen file and read meta data
fid = fopen(filename);
fgets(fid);
buff = fgets(fid);
N = sscanf(buff,'%*s %*s %*s %i %*s %*s %*s %i %*s %*s %*s %i');
Nk=N(1); Nb=N(2); Ni=N(3);
if have_charge
    Nentries = Ni+1;
else
    Nentries = Ni;
end

%% read data
k = zeros(3,Nk);
w = zeros(1,Nk);

for ck=1:Nk
    % read k-point and weight
    buff=[];
    while isempty(strfind(buff,'k-point'))
        buff=fgets(fid);
        if feof(fid)
            error('unexpected end of file')
        end
    end
    dat = sscanf(buff,'%*s %*i %*s %lf %lf %lf %*s %*s %lf');
    k(:,ck) = dat(1:3);
    w(ck) = dat(4);
    
    % initialize data
    B(ck).E = zeros(1,Nb);
    B(ck).occ = zeros(1,Nb);
    B(ck).s = zeros(Nentries,Nb);
    B(ck).p = zeros(Nentries,Nb);
    B(ck).d = zeros(Nentries,Nb);
    B(ck).tot = zeros(Nentries,Nb);
    if have_f
        B(ck).f = zeros(Nentries,Nb);
    end
    if lm_decomp
        B(ck).px = zeros(Nentries,Nb);
        B(ck).py = zeros(Nentries,Nb);
        B(ck).pz = zeros(Nentries,Nb);
        B(ck).dxy = zeros(Nentries,Nb);
        B(ck).dyz = zeros(Nentries,Nb);
        B(ck).dxz = zeros(Nentries,Nb);
        B(ck).dz2 = zeros(Nentries,Nb);
        B(ck).dx2 = zeros(Nentries,Nb);
        if have_f
            B(ck).fy3x2 = zeros(Nentries,Nb);
            B(ck).fxyz = zeros(Nentries,Nb);
            B(ck).fyz2 = zeros(Nentries,Nb);
            B(ck).fz3 = zeros(Nentries,Nb);
            B(ck).fxz2 = zeros(Nentries,Nb);
            B(ck).fzx2 = zeros(Nentries,Nb);
            B(ck).fx3 = zeros(Nentries,Nb);
        end
    end

    for cb = 1:Nb
        % read energy and occupation
        buff=[];
        while isempty(strfind(buff,'band'))
            buff=fgets(fid);
            if feof(fid)
                error('unexpected end of file')
            end
        end
        dat = sscanf(buff,'%*s %*i %*s %*s %lf %*s %*s %lf');
        B(ck).E(cb) = dat(1);
        B(ck).occ(cb) = dat(2);

        % read projections for each ion
        buff=[];
        while isempty(strfind(buff,'ion'))
            buff=fgets(fid);
            if feof(fid)
                error('unexpected end of file')
            end
        end
        for ci=1:Ni
            buff = fgets(fid);
 
            dat = sscanf(buff,['%*i ' format]);
            
            B(ck).s(ci,cb) = dat(1);
            B(ck).tot(ci,cb) = dat(end);
            if ~lm_decomp
                B(ck).p(ci,cb) = dat(2);
                B(ck).d(ci,cb) = dat(3);
                if have_f
                    B(ck).f(ci,cb) = dat(4);
                end
            else
                B(ck).py(ci,cb) = dat(2);
                B(ck).pz(ci,cb) = dat(3);
                B(ck).px(ci,cb) = dat(4);
                B(ck).p(ci,cb) = sum(dat(2:4));
                B(ck).dxy(ci,cb) = dat(5);
                B(ck).dyz(ci,cb) = dat(6);
                B(ck).dxz(ci,cb) = dat(7);
                B(ck).dz2(ci,cb) = dat(8);
                B(ck).dx2(ci,cb) = dat(9);
                B(ck).d(ci,cb) = sum(dat(5:9));
                if have_f
                    B(ck).fy3x2(ci,cb) = dat(10);
                    B(ck).fxyz(ci,cb) = dat(11);
                    B(ck).fyz2(ci,cb) = dat(12);
                    B(ck).fz3(ci,cb) = dat(13);
                    B(ck).fxz2(ci,cb) = dat(14);
                    B(ck).fzx2(ci,cb) = dat(15);
                    B(ck).fx3(ci,cb) = dat(16);
                    B(ck).f(ci,cb) = sum(dat(10:16));
                end
            end
        end
        
        % read charge, LORBIT = 12 only
        if have_charge
            buff=[];
            while isempty(strfind(buff,'charge'))
                buff=fgets(fid);
                if feof(fid)
                    error('unexpected end of file')
                end
            end
            dat = sscanf(buff,['%*s ' format]);
            B(ck).s(end,cb) = dat(1);
            B(ck).tot(end,cb) = dat(end);
            B(ck).py(end,cb) = dat(2);
            B(ck).pz(end,cb) = dat(3);
            B(ck).px(end,cb) = dat(4);
            B(ck).p(end,cb) = sum(dat(2:4));
            B(ck).dxy(end,cb) = dat(5);
            B(ck).dyz(end,cb) = dat(6);
            B(ck).dz2(end,cb) = dat(7);
            B(ck).dxz(end,cb) = dat(8);
            B(ck).dx2(end,cb) = dat(9);
            B(ck).d(end,cb) = sum(dat(5:9));
            if have_f
                B(ck).fy3x2(end,cb) = dat(10);
                B(ck).fxyz(end,cb) = dat(11);
                B(ck).fyz2(end,cb) = dat(12);
                B(ck).fz3(end,cb) = dat(13);
                B(ck).fxz2(end,cb) = dat(14);
                B(ck).fzx2(end,cb) = dat(15);
                B(ck).fx3(end,cb) = dat(16);
                B(ck).f(end,cb) = sum(dat(10:16));
            end
            
        end
    end
end

fclose(fid);

end
