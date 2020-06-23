function M = read_bin_mat(filename)
% function to read matrix from binary file
% CALL: M = READ_BIN_MAT(filename)
%
% input args
% - filename: name of the binary file
%
% output args
% - M: output matrix

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


fid = fopen(filename);
head = fread(fid,5,'char');
dim = fread(fid,1,'uint64');
S = fread(fid,dim,'uint64')';

if (all([102 77 97 116 0]'==head))
    % fMat
    M = fread(fid,prod(S),'double');
end
if (all([99 77 97 116 0]'==head))
    % cMat
    blob = fread(fid,2*prod(S),'double');
    M = blob(1:2:2*prod(S)-1) + 1i*blob(2:2:2*prod(S));
end

% reshape according to S
M = reshape(M,S);

fclose(fid);
