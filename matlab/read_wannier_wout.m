function [B,R,Wp,s,Apf,Apc,T,id,Na,Nw] = read_wannier_wout(filename)
% function to read data from wannier90 output
% CALL: [B,R,Wp,s,Apf,Apc,T,id,Na,Nw] = READ_WANNIER_WOUT('wannier90.wout')
%
% input args:
% - filename: name of the wannier90 wout file, default: 'wannier90.wout'
%
% output args:
% - B: lattice vectors as columns, size 3x3
% - R: reciprocal space vectors as columns, size 3x3
% - Wp: wannier centers as columns, size 3xNw
% - s: spread for each wannier function, size 1xNw
% - Apf: atomic positions in fractional as columns, size 3xNa
% - Apc: atomic positions in cartesian as columns, size 3xNa
% - T: atomic type, size 1xNa
% - id: atomic identifiers, size 1xunique(Na)
% - Na: # atomic positions
% - Nw: # wannier functions

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<1
    filename = 'wannier90.wout';
end

%% open file
fid = fopen(filename);

%% read B, R
buff=[];
while isempty(strfind(buff,'Lattice Vectors'))
    buff = fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
B = zeros(3);
for cb=1:3
    buff = fgets(fid);
    B(:,cb) = sscanf(buff,'%*s %lf %lf %lf');
end

buff=[];
while isempty(strfind(buff,'Reciprocal-Space Vectors'))
    buff = fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
R = zeros(3);
for cr=1:3
    buff = fgets(fid);
    R(:,cr) = sscanf(buff,'%*s %lf %lf %lf');
end

%% read Apf, Apc, T, id, Na
buff=[];
while isempty(strfind(buff,'Fractional Coordinate'))
    buff = fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
fgets(fid);

Apf=[]; Apc=[]; T = []; id = {};
buff = fgets(fid);
c=1;
while isempty(strfind(buff,'*----'))
    if feof(fid)
        error('unexpected end of file')
    end
    id{c} = sscanf(buff,'%*s %s',1);
    c=c+1;
    
    % Apf Apc
    Apf = [Apf sscanf(buff,'%*s %*s %*i %lf %lf %lf',3)];
    Apc = [Apc sscanf(buff,'%*s %*s %*i %*lf %*lf %*lf %*s %lf %lf %lf',3)];
    
    buff = fgets(fid);
end

Na = size(Apf,2);

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

%% read Wp, s, Nw
buff=[];
while isempty(strfind(buff,'Number of Wannier Functions'))
    buff = fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end
Nw = sscanf(buff,'%*s %*s %*s %*s %*s %*s %lf',1);


buff=[];
while isempty(strfind(buff,'Final State'))
    buff = fgets(fid);
    if feof(fid)
        error('unexpected end of file')
    end
end

Wp = zeros(3,Nw); s=zeros(1,Nw);
for cw=1:Nw
    buff = fgets(fid);
    lfbuff = sscanf(buff,'%*s %*s %*s %*s %*lf %*c %lf %*c %lf %*c %lf %*c %lf',4);
    Wp(:,cw) = lfbuff(1:3)';
    s(cw) = lfbuff(4);
end

fclose(fid);

end
