function [H,R] = read_bin_sparse_hr(filename)
% function to read binary sparse hr from file
% CALL: [H,R] = READ_BIN_SPARSE_HR(filename)
%
% input args
% - filename: name of the binary file
%
% output args
% - H: cell array of size 1xNR
% - R: matrix of R vectors, size 3xNR

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


% open file
fid = fopen(filename);

% read meta data
head = fread(fid,5,'char');
dim = fread(fid,1,'uint32');
Nw = fread(fid,1,'uint32');
NR = fread(fid,1,'uint32');

R = zeros(dim,NR); H = cell(1,NR);

% single precision, short, real
if strcmp('hrssr',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'float64');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint16');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            h(c) = fread(fid,1,'float32');
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% single precision, long, real
if strcmp('hrlsr',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'float64');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint32');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            h(c) = fread(fid,1,'float32');
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% single precision, short, complex
if strcmp('hrssc',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'double');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint16');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            tmp = fread(fid,2,'float32');
            h(c) = tmp(1) + 1i*tmp(2);
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% single precision, long, complex
if strcmp('hrlsc',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'double');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint32');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            tmp = fread(fid,2,'float32');
            h(c) = tmp(1) + 1i*tmp(2);
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% double precision, short, real
if strcmp('hrsdr',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'float64');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint16');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            h(c) = fread(fid,1,'float64');
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% double precision, long, real
if strcmp('hrldr',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'float64');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint32');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            h(c) = fread(fid,1,'float64');
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% double precision, short, complex
if strcmp('hrsdc',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'double');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint16');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            tmp = fread(fid,2,'float64');
            h(c) = tmp(1) + 1i*tmp(2);
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end

% double precision, long, complex
if strcmp('hrldc',char(head'))
    for cr=1:NR
        R(:,cr) = fread(fid,3,'double');
        N = fread(fid,1,'uint64');
        
        m = zeros(1,N); n = zeros(1,N); h = zeros(1,N);
        for c=1:N
            ci = fread(fid,2,'uint32');
            m(c) = ci(1)+1; n(c) = ci(2)+1;
            tmp = fread(fid,2,'float64');
            h(c) = tmp(1) + 1i*tmp(2);
        end
        H{cr} = sparse(m,n,h,Nw,Nw);
    end
    fclose(fid);
    return
end


fclose(fid);

% complain
error(['bad header: ''' head ''''])

end
