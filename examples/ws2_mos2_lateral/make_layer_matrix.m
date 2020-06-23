%% read unit cell from wbh, get positions in cartesian and shift to origin
[B,Ap,id] = cell_from_wbh('wannier90.wbh');
Ap = B*Ap;
Ap(2,:) = Ap(2,:) - min(Ap(2,:));

%% get indices for each type
Im = []; Iw = []; Is = [];
for c=1:numel(id)
    if ~isempty(strfind(id{c},'Mo'))
        Im = [Im c];
    end
    if ~isempty(strfind(id{c},'W'))
        Iw = [Iw c];
    end
    if ~isempty(strfind(id{c},'S'))
        Is = [Is c];
    end
end

%% plot and data cursor
% close all
% hold on
% scatter3(Ap(1,Iw),Ap(2,Iw),Ap(3,Iw),'blue'), axis equal
% scatter3(Ap(1,Im),Ap(2,Im),Ap(3,Im),'green'), axis equal
% scatter3(Ap(1,Is),Ap(2,Is),Ap(3,Is),'black'), axis equal


%% midsection should include 6 Mo atoms and 6 W atoms
ub = 41.5; lb = 10;
Jmid = Ap(1,:)<ub & Ap(1,:)>lb;

%% repeat first 2 W, last 2 Mo
ub = 13.82; lb = 10;
Jw = Ap(1,:)<ub & Ap(1,:)>lb;
ub = 41.45; lb = 37.75;
Jm = Ap(1,:)<ub & Ap(1,:)>lb;

%% append midsection, repeat N times on each side
N = 33; % results in ~400A device length
sh = Ap(1,Iw(2)) - Ap(1,Iw(1));

Wtail = []; Wtail_id = {};
for c=1:N
    Wtail = [Ap(:,Jw)-[1;0;0]*ones(1,sum(Jw))*c*sh Wtail];
    Wtail_id = [id(Jw) Wtail_id];
end
Mtail = []; Mtail_id = {};
for c=1:N
    Mtail = [Mtail Ap(:,Jm)+[1;0;0]*ones(1,sum(Jm))*c*sh];
    Mtail_id = [Mtail_id id(Jm)];
end

%% final structure is tails on each side and midsection
fAp = [Wtail Ap(:,Jmid) Mtail];
fAp(1,:) = fAp(1,:) - min(fAp(1,:));
fid = [Wtail_id id(Jmid) Mtail_id];

%% duplicate 5 times along z to have nn only interactions
zh = B(3,1);
fAp0 = fAp; fAp0(3,:) = fAp0(3,:) + 0*zh;
fAp1 = fAp; fAp1(3,:) = fAp1(3,:) + 1*zh;
fAp2 = fAp; fAp2(3,:) = fAp2(3,:) + 2*zh;
fAp3 = fAp; fAp3(3,:) = fAp3(3,:) + 3*zh;
fAp4 = fAp; fAp4(3,:) = fAp4(3,:) + 4*zh;

fAp = [fAp0 fAp1 fAp2 fAp3 fAp4];
fid = [fid fid fid fid fid];

[~,I] = sortrows(fAp');
fAp = fAp(:,I);
fid = fid(I);

fB = diag([Ap(1,45)-Ap(1,33)+2*N*sh max(Ap(2,:))+3 5*B(3,1)]);


write_lattice_dat(fB,fAp,fid,6,3.5,'Layer_Matrix.dat')


