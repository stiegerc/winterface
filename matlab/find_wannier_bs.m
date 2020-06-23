function [E,k,kpl,N,ticks] = find_wannier_bs(H,R,B,re,kpt,N,tol,bL)
% function to compute bandstructure from wannier90 hr data
% CALL: [E,k,kpl,N,ticks] = FIND_WANNIER_BS(H,R,B,re,kpt,N,tol,bL)
%
% input args:
% - H: hamiltonian cellarray as produced by 'read_wannier_hr'
% - R: R vectors matching R as produced by 'read_wannier_hr'
% - B: base matching hamiltonian data as produced by 'read_wannier_wout'
% - re: bool to discard imaginary part of hamiltonian data, default: false
% - kpt: points in reciprocal k space to produce k paths, size 3xNkpt,
%           default [-.5 0 0; .5 0 0]'
% - N: # k points along each line, size 1 for autoscaling or 1x(Nkpt-1),
%           default: 300
% - tol: tolerance below which data in H is considered 0, default: 0
% - bL: list of blacklisted wannier functions, size 1x(<Nw), default []
%
% output args:
% - E: bandstructure data, size NwxNk
% - k: k points, size 3xNk
% - kpl: k plotting vector, size NwxNk
% - N: number of k points per trace, size 1x(Nkpt-1)
% - ticks: ticks according to kpt, size 1xNkpt

% 2014-2019, ETH Zurich, Integrated Systems Laboratory
% Authors: Christian Stieger


%% default argument
if nargin<4
    re = false;
end
if nargin<5
    kpt = [-.5 0 0; .5 0 0]';
end
if nargin<6
    N = 300;
end
if nargin<7
    tol = 0;
end
if nargin<8
    bL=[];
end


%% bitching
if size(kpt,1)~=3
    error('need size(kpt,1)==3')
end


%% treat H according to re, tol, bL
Nw = size(H{1},1); NR = size(R,2);
wLH = true(1,NR); wLw = true(1,Nw); wLw(bL) = false;
for cR=1:NR
    if re
        H{cR} = real(H{cR});
    end
    H{cR} = H{cR}(wLw,wLw);
    
    I = abs(H{cR})>=tol;
    if any(any(I))
        H{cR}(~I) = 0;
    else
        wLH(cR) = false;
    end
end
R = R(:,wLH);
H = H(wLH);


%% get N, kpl, ticks
Nw = size(H{1},1); NR = size(R,2);

l = zeros(1,size(kpt,2)-1);
Bhat = 2*pi*inv(B)';
for c=1:length(l)
    l(c) = norm(Bhat*(kpt(:,c+1)-kpt(:,c)));
end
if length(N)~=length(l)
    N = ceil(l/sum(l)*N);
end

kpl = []; ticks = 0;
for c=1:length(l)
    kpl = [kpl linspace(ticks(end),ticks(end)+l(c),N(c))];
    ticks = [ticks ticks(end)+l(c)];
end
kpl = ones(Nw,1)*kpl;

%% get k, Nk
k = kpath(kpt,N);
Nk = size(k,2);


%% get E
E = zeros(Nw,Nk);
parfor ck=1:Nk
    Ek = zeros(Nw);
    for cR=1:NR
        Ek = Ek + H{cR}*exp(-1i*2*pi*R(:,cR)'*k(:,ck));
    end
    
    E(:,ck) = sort(real(eig(Ek)));
end

end
