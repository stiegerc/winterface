% plot1: 3D plot of the bandstructure, including the trace of plot2 in red
% plot2: bandstructure along the Gamma-K-M-Gamma trace,
%    colored lines: exact reference
%    dashed black lines: approximate bandstructure using a cutoff value
% plot3: 3D plot of the atomic unit cell, including wannier centers in red
%         Mo: blue, S: green

fs = 20; % font size

close all
figure

prefix = 'output/w/';
vb_I = 1:7;
cb_I = 8:11;

I = [vb_I(end-3:end) cb_I(1:4)];
Ef = read_outcar('data/OUTCAR');
[~,R] = read_wannier_wout('data/wannier90.wout');
kpts = [
    0.00000 0.00000 0.00000 % Gamma
    0.33333 0.33333 0.00000 % K
    0.50000 0.00000 0.00000 % M
    0.00000 0.00000 0.00000 % Gamma
]';
lbls = {'\Gamma','K','M','\Gamma'};

Emin = -4; Emax = 4.7;


%% plot mesh bs
Emesh_sc = read_bin_mat([prefix 'Escal_mesh.bin']);

[~,Ny,Nx] = size(Emesh_sc);
[kx,ky] = meshgrid(linspace(-.5,.5,Nx),linspace(-.5,.5,Ny));
for c=1:numel(kx)
    % convert mesh to cartesian
    cc = R * [kx(c) ; ky(c) ; 0];
    kx(c) = cc(1);
    ky(c) = cc(2);
end

kpts_ = R*kpts;
lx = [kpts_(1,1) kpts_(1,1) kpts_(1,2) kpts_(1,2) kpts_(1,3) kpts_(1,3)];
ly = [kpts_(2,1) kpts_(2,1) kpts_(2,2) kpts_(2,2) kpts_(2,3) kpts_(2,3)];
lz = [Emin Emax Emin Emax Emin Emax];

subplot(1,3,1), title('3D bandstructure')
hold on
for cb=I
    e = squeeze(Emesh_sc(cb,:,:))-Ef;
    surf(kx,ky,e)
end
plot3(lx(1:2),ly(1:2),lz(1:2),'red','LineWidth',2)
plot3(lx(3:4),ly(3:4),lz(3:4),'red','LineWidth',2)
plot3(lx(5:6),ly(5:6),lz(5:6),'red','LineWidth',2)
scatter3(lx(1:2),ly(1:2),lz(1:2),'filled','red')
scatter3(lx(3:4),ly(3:4),lz(3:4),'filled','red')
scatter3(lx(5:6),ly(5:6),lz(5:6),'filled','red')
plot3([lx(1) lx(3)],[ly(1) ly(3)],[Emin Emin],'red','LineWidth',2)
plot3([lx(3) lx(5)],[ly(3) ly(5)],[Emin Emin],'red','LineWidth',2)
plot3([lx(5) lx(1)],[ly(5) ly(1)],[Emin Emin],'red','LineWidth',2)
plot3([lx(1) lx(3)],[ly(1) ly(3)],[Emax Emax],'red','LineWidth',2)
plot3([lx(3) lx(5)],[ly(3) ly(5)],[Emax Emax],'red','LineWidth',2)
plot3([lx(5) lx(1)],[ly(5) ly(1)],[Emax Emax],'red','LineWidth',2)
hold off
xlabel('k_x'), ylabel('k_y'), zlabel('E-E_F [eV]')
box on
grid on
axis([-2 2 -2 2 Emin Emax])
shading faceted
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

delta = 0.25;
for c=1:3
    t = text(kpts_(1,c),kpts_(2,c),Emax+delta,lbls{c});
    set(t,'FontSize',fs)
    set(t,'Color','red')
    set(t,'FontWeight','bold')
end
for c=1:3
    t = text(kpts_(1,c),kpts_(2,c),Emin-delta,lbls{c});
    set(t,'FontSize',fs)
    set(t,'Color','red')
    set(t,'FontWeight','bold')
end



%% plot trace bs
Efold = read_bin_mat([prefix 'Efold.bin']);
Escal = read_bin_mat([prefix 'Escal.bin']);
pos = read_bin_mat([prefix 'pos.bin']);
path = read_bin_mat([prefix 'path.bin']);

[Nb,Nk] = size(Escal);

% find indices for labels
J = zeros(1,4);
for c=2:3
    n = zeros(1,Nk);
    for d=1:Nk
        n(d) = norm(path(:,d)-kpts(:,c));
    end
    [~,j] = min(n);
    J(c) = j;
end
J(1) = 1; J(end) = Nk;

subplot(1,3,2), title('bandstructure trace')
hold on
for cb=I
    plot(pos,Escal(cb,:)-Ef,'LineWidth',3)
    plot(pos,Efold(cb,:)-Ef,'k--','LineWidth',1)
end
hold off
set(gca,'XTick',pos(J));
set(gca,'XTickLabel',lbls);
box on
grid on
ylabel('E-E_F [eV]')
axis([pos(1) pos(end), Emin Emax])
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)



%% plot unit cell
[B,Ap1] = read_poscar([prefix 'iwp_0001_[Mo]_05_(0.063,0.113).psc']);
[~,Ap2] = read_poscar([prefix 'iwp_0002_[S]_03_(0.061,0.099).psc']);
[~,Ap3] = read_poscar([prefix 'iwp_0003_[S]_03_(0.061,0.099).psc']);
Ap1 = B*Ap1;
Ap2 = B*Ap2;
Ap3 = B*Ap3;
Ap1(3,:) = Ap1(3,:)-Ap2(3,1);
Ap3(3,:) = Ap3(3,:)-Ap2(3,1);
Ap2(3,:) = Ap2(3,:)-Ap2(3,1);

B(3,3) = 3.5;
P1 = [
    0 0 0
    1 0 0
    1 1 0
    0 1 0
    0 0 0
]';
P2 = [
    0 0 1
    1 0 1
    1 1 1
    0 1 1
    0 0 1
]';
P3 = [
    0 0 0
    0 0 1
    1 0 1
    1 0 0
    1 1 0
    1 1 1
    0 1 1
    0 1 0
]';

subplot(1,3,3), title('unit cell')
hold on
scatter3(Ap1(1,1),Ap1(2,1),Ap1(3,1),'filled','blue','SizeData',120)
scatter3(Ap1(1,2:end),Ap1(2,2:end),Ap1(3,2:end),'filled','red','SizeData',40)
scatter3(Ap2(1,1),Ap2(2,1),Ap2(3,1),'filled','green','SizeData',120)
scatter3(Ap2(1,2:end),Ap2(2,2:end),Ap2(3,2:end),'filled','red','SizeData',40)
scatter3(Ap3(1,1),Ap3(2,1),Ap3(3,1),'filled','green','SizeData',120)
scatter3(Ap3(1,2:end),Ap3(2,2:end),Ap3(3,2:end),'filled','red','SizeData',40)
plot3([Ap3(1,1) Ap1(1,1) Ap2(1,1)],[Ap3(2,1) Ap1(2,1) Ap2(2,1)],[Ap3(3,1) Ap1(3,1) Ap2(3,1)],'green','LineWidth',2)
pp = B*P1; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P2; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P3; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
hold off
axis equal
box on
grid on
xlabel('x'), ylabel('y'), zlabel('z')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)
