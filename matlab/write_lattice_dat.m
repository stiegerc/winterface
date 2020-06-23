function write_lattice_dat(B,Ap,id,nn,bl,filename)
% function to write OMEN lattice_dat file
% CALL: WRITE_LATTICE_DAT(B,Ap,id,nn,bl,filename)
%
% input args:
% - B: basis vectors for the lattice
% - Ap: atomic positions
% - id: cell array of identifier strings for each position in Ap
% - nn: number of next neighbours
% - bl: bond length
% - filename: a name for the output file

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% bitching
if size(B,1)~=size(B,2)
    error('B must be square')
end
if size(B,1)~=size(Ap,1)
    error('B and Ap dimension mismatch')
end
if numel(id)~=size(Ap,2)
    error('need size(Ap.2)==numel(id)')
end

%% function body
fid = fopen(filename,'w');
fprintf(fid,'%i %i %i %i %i \n\n',[size(Ap,2) nn 0 0 0]);
fprintf(fid,'%16.12f %16.12f %16.12f \n',B(:,1));
fprintf(fid,'%16.12f %16.12f %16.12f \n',B(:,2));
fprintf(fid,'%16.12f %16.12f %16.12f \n\n',B(:,3));

for c=1:size(Ap,2)
    fprintf(fid,'%5s %16.12f %16.12f %16.12f \n',id{c},Ap(:,c));
end

fprintf(fid,'\n%f',bl);

fclose(fid);

end
