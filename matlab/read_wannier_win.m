function [B,Ap,T,id,mpg,kpt,Na,Nw,Nb] = read_wannier_win(filename)
% function to read data from wannier90 input
% CALL: [B,Ap,T,id,mpg,kpt,Na,Nw,Nb] = READ_WANNIER_WIN('wannier90.win')
%
% input args:
% - filename: name of the wannier90 win file, default: 'wannier90.win'
%
% output args:
% - B: lattice vectors as columns, size 3x3
% - Ap: atomic positions in cartesian, size 3xNa
% - T: atomic type, size 1xNa
% - id: atomic identifiers, size 1xunique(Na)
% - mpg: mp_grid the kpoint mesh parameters, size 1x3
% - kpt: kpoints, size 3xNk
% - Na: # atomic positions
% - Nw: # wannier functions
% - Nb: # bands

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'wannier90.win';
end

%% open file
fid = fopen(filename);

%% set output to empty initially
B=[]; Ap=[]; T=[]; id=[]; mpg=[]; kpt=[]; Na=[]; Nw=[]; Nb=[];

%% parse file
while ~feof(fid)
    buff = fgets(fid);
    if buff(1)=='#' || buff(1)=='!'
        continue
    end
    
    if ~isempty(strfind(buff,'num_bands'))
        Nb = sscanf(buff,'%*s %*c %i',1);
    end
    if ~isempty(strfind(buff,'num_wann'))
        Nw = sscanf(buff,'%*s %*c %i',1);
    end
    
    if ~isempty(strfind(buff,'begin unit_cell_cart'))
        B = zeros(3,3);
        B(:,1) = sscanf(fgets(fid),'%lf %lf %lf',3);
        B(:,2) = sscanf(fgets(fid),'%lf %lf %lf',3);
        B(:,3) = sscanf(fgets(fid),'%lf %lf %lf',3);
    end
    
    if ~isempty(strfind(buff,'begin atoms_cart'))
        Ap = [];
        id = {};
        
        % read Ap, all ids
        buff = fgets(fid);
        while isempty(strfind(buff,'end atoms_cart'))
            Ap = [Ap sscanf(buff,'%*s %lf %lf %lf',3)];
            id{end+1} = sscanf(buff,'%s',1);
            buff = fgets(fid);
        end

        Na = size(Ap,2);
        
        % find T, unique ids
        uid = unique(id);
        T = zeros(1,Na);
        for c1=1:length(id)
            for c2=1:length(uid)
                if strcmp(id{c1},uid{c2})
                    T(c1)=c2;
                end
            end
        end
        id=uid;
    end
    
    if ~isempty(strfind(buff,'mp_grid'))
        mpg = sscanf(buff,'%*s %*c %i %i %i')';
    end
    
    if ~isempty(strfind(buff,'begin kpoints'))
        kpt = [];
        buff = [];
        while isempty(strfind(buff,'end kpoints'))
            buff = fgets(fid);
            kpt = [kpt sscanf(buff,'%lf %lf %lf',3)];
        end
    end
end

fclose(fid);

end
