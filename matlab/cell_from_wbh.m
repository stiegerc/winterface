function [B,Ap,id] = cell_from_wbh(filename)
% function to read unit cell from a wbh
% CALL: [B,Ap,id] = CELL_FROM_WBH('wbh.wad')
%
% terminology:
% - Na: # of atoms in the unit cell
%
% input args:
% - filename: name of the wbh file, default: 'wbh.wad'
%
% output args:
% - B: cell vectors as columns, size 3x3
% - Ap: atomic positions in basis B, size 3xNa
% - id: cell array of strings, size 1xNa

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default filename
if nargin<1
    filename = 'wbh.wad';
end

%% open file, check header
fid = fopen(filename);
head = fread(fid,5,'char');
if (~all([119 97 100 57 48]'==head))
    disp('bad header')
    return
end

%% read cell dimensions
D = fread(fid,1,'uint32');
N = fread(fid,1,'uint32');

%% read B and Ap
B = zeros(D,D);
B(:) = fread(fid,D*D,'double');
Ap = zeros(D,N);
Ap(:) = fread(fid,D*N,'double');

%% read id
id = cell(1,N);
for c=1:N
    id{c} = fread(fid,fread(fid,1,'uint8'),'char');
    id{c} = char(id{c}(1:end-1)');
end

end
