% wannier bonds hamiltonian class

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


classdef wbh
    properties
        cell
        Norb
        I
        Ef
    end
    methods
        %% ctor from file
        function obj = wbh(filename)
            % default filename
            if nargin<1
                filename = 'wannier90.wbh';
            end
            
            % open file, check header
            fid = fopen(filename);
            head = fread(fid,5,'char');
            if (~all([119 97 100 57 48]'==head))
                disp('bad header')
                return
            end
            
            % read cell dimensions
            D = fread(fid,1,'uint32');
            N = fread(fid,1,'uint32');
            
            % read B and Ap
            obj.cell.B = zeros(D,D);
            obj.cell.B(:) = fread(fid,D*D,'double');
            obj.cell.Ap = zeros(D,N);
            obj.cell.Ap(:) = fread(fid,D*N,'double');
            
            % read id
            obj.cell.id = cell(1,N);
            for c=1:N
                obj.cell.id{c} = fread(fid,fread(fid,1,'uint8'),'char');
                obj.cell.id{c} = char(obj.cell.id{c}(1:end-1)');
            end
            
            % read Norb
            obj.Norb = fread(fid,N,'uint32');
            
            % read Nb
            Nb = fread(fid,1,'uint32');
            
            % read index pairs, R vectors, H matrices
            for c=1:Nb
                dat = fread(fid,3,'uint32');
                obj.I(c).i1 = dat(1)+1;
                obj.I(c).i2 = dat(2)+1;
                NR = dat(3);
                
                obj.I(c).R = zeros(D,NR);
                obj.I(c).R(:) = fread(fid,D*NR,'double');
                
                obj.I(c).H = cell(1,NR);
                for d=1:NR
                    obj.I(c).H{d} = zeros(obj.Norb(obj.I(c).i1), ...
                                          obj.Norb(obj.I(c).i2));
                    blob = fread(fid,2*numel(obj.I(c).H{d}),'double');
                    obj.I(c).H{d}(:) = blob(1:2:end-1) + 1i*blob(2:2:end);
                end
            end
            
            % read Ef
            obj.Ef = fread(fid,1,'double');
            
            % close file
            fclose(fid);
        end        
            
        %% information
        function j = cbegin(obj,i1)
            for c=1:numel(obj.I)
                if (obj.I(c).i1 == i1)
                    j = c;
                    return
                end
            end
            j = numel(obj.I);
        end
        function j = cend(obj,i1)
            for c=cbegin(obj,i1):numel(obj.I)
                if (obj.I(c).i1 ~= i1)
                    j = c-1;
                    return
                end
            end
            j = numel(obj.I);
        end
        function j = ind(obj,i1,i2)
            for c=cbegin(obj,i1):cend(obj,i1)
                if (obj.I(c).i2 == i2)
                    j = c;
                    return
                end
            end
            j = NaN;
        end
        function pb = pbond(obj,i1,i2)
            pb = obj.cell.Ap(:,i2)-obj.cell.Ap(:,i1);
        end
        function RR = range(obj)
            RR = zeros(1,3);
            for c=1:numel(obj.I)
                tmp = [RR ; obj.I(c).R'];
                RR = unique(sortrows(tmp),'rows');
            end
            RR = RR';
        end
        function T = id_to_T(obj,id)
            T = -ones(1,numel(id));
            for c=1:numel(id)
                for d=1:numel(obj.cell.id)
                    if strcmp(obj.cell.id{d},id{c})
                        T(c)=d;
                    end
                end
            end
        end
        
        %% get interactions
        function [h,i,j] = dget_h(obj,i1,i2,b,tol)
            i = ind(obj,i1,i2);
            if isnan(i)
                h = [];
                return
            end
            b = b - obj.pbond(i1,i2);
            
            NR = size(obj.I(i).R,2);
            for c=1:NR
                if all(abs(obj.I(i).R(:,c) - b)<tol)
                    j = c;
                    h = obj.I(i).H{j};
                    return
                end
            end
            
            j = NaN;
            h = [];
            
        end
        function S = smat(obj,I1,I2,R)
            S = zeros(numel(I1),numel(I2));
            for c1=1:numel(I1)
                for c2=1:numel(I2)
                    i1 = I1(c1); i2 = I2(c2); pb = obj.pbond(i1,i2);
                    h = obj.dget_h(i1,i2,pb+R,1e-3);
                    S(c1,c2) = sum(svd(h));
                end
            end
        end
        
        %% write to file
        function write(obj,filename)
            % open file
            fid = fopen(filename,'w');
            
            % write header
            fwrite(fid,[119 97 100 57 48],'char');
            
            % write cell
            fwrite(fid,size(obj.cell.B,1),'uint32');
            fwrite(fid,size(obj.cell.Ap,2),'uint32');
            fwrite(fid,obj.cell.B,'double');
            fwrite(fid,obj.cell.Ap,'double');
            for c=1:numel(obj.cell.id)
                fwrite(fid,numel(obj.cell.id{c})+1,'uint8');
                fwrite(fid,obj.cell.id{c},'char');
                fwrite(fid,0,'char');
            end
            
            % write Norb
            fwrite(fid,obj.Norb,'uint32');
            
            % write #pairs
            fwrite(fid,numel(obj.I),'uint32');

            % write #pairs, pairs, #bonds, R and H
            for c=1:numel(obj.I)
                fwrite(fid,[obj.I(c).i1-1 obj.I(c).i2-1 numel(obj.I(c).H)],'uint32');
                fwrite(fid,obj.I(c).R,'double');
                for d=1:numel(obj.I(c).H)
                    re = real(obj.I(c).H{d}); im = imag(obj.I(c).H{d}); 
                    buff = zeros(2*numel(re),1);
                    buff(1:2:end-1) = re(:); buff(2:2:end) = im(:);
                    fwrite(fid,buff,'double');
                end
            end
            
            % write Ef
            fwrite(fid,obj.Ef,'double');
            
            % close file
            fclose(fid);
        end
        function write_poscar(obj,filename)
            write_poscar(obj.cell.B,obj.cell.Ap,1,ones(1,size(obj.cell.Ap,2)),...
                obj.cell.id,true,filename)
        end
    end
end
