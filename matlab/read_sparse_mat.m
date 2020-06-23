function H = read_sparse_mat(filename)
% function to read sparse matrix from file, omen style
% CALL: H = READ_SPARSE_MAT(filename)
%
% input args:
% - filename: name of the sparse matrix file
%
% output args:
% - H: sparse matrix

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


fid = fopen(filename);
head = fread(fid,3,'double');

blob = fread(fid,4*head(2),'double');
blob = reshape(blob,4,head(2))';
fclose(fid);


if head(3)==0
    H = sparse(blob(:,1)+1,blob(:,2)+1,blob(:,3)+1i*blob(:,4),head(1),head(1));
else
    H = sparse(blob(:,1),blob(:,2),blob(:,3)+1i*blob(:,4),head(1),head(1));
end

end
