% figure 1
% plot1: 3D plot of the bandstructure, including the trace of plot2 in red
% plot2: bandstructure along the Gamma-X-S-Y-Gamma trace
%
% figure 2
% plot1: bandstructure along the Gamma-X-S-Y-Gamma trace
% plot2: first derivative of plot1 in direction of the trace
% plot3: second derivative of plot1 in direction of the trace


fs = 20; % font size

close all
figure

prefix = 'output/b/';
vb_I = 1:14;
cb_I = 15:22;

I = [vb_I(end-7:end) cb_I(1:8)];
Ef = read_outcar('data/OUTCAR');
B = dlmread([prefix 'NB.mat']); R = 2*pi*inv(B');
kpts = [
    0.00000 0.00000 0.00000 % Gamma
    0.50000 0.00000 0.00000 % X
    0.50000 0.50000 0.00000 % S
    0.00000 0.50000 0.00000 % Y
    0.00000 0.00000 0.00000 % Gamma
]';
lbls = {'\Gamma','X','S','Y','\Gamma'};

Emin = -4; Emax = 4.7;


%% fig2, plot mesh bs
Emesh = read_bin_mat([prefix 'Efold_mesh.bin']);

[~,Ny,Nx] = size(Emesh);
[kx,ky] = meshgrid(linspace(-.5,.5,Nx),linspace(-.5,.5,Ny));
for c=1:numel(kx)
    % convert mesh to cartesian
    cc = R * [kx(c) ; ky(c) ; 0];
    kx(c) = cc(1);
    ky(c) = cc(2);
end

kpts_ = R*kpts;
lx = [kpts_(1,1) kpts_(1,1) kpts_(1,2) kpts_(1,2) kpts_(1,3) kpts_(1,3) kpts_(1,4) kpts_(1,4)];
ly = [kpts_(2,1) kpts_(2,1) kpts_(2,2) kpts_(2,2) kpts_(2,3) kpts_(2,3) kpts_(2,4) kpts_(2,4)];
lz = [Emin Emax Emin Emax Emin Emax Emin Emax];

subplot(1,2,1), title('3D bandstructure')
hold on
for cb=I
    e = squeeze(Emesh(cb,:,:))-Ef;
    surf(kx,ky,e)
end
plot3(lx(1:2),ly(1:2),lz(1:2),'red','LineWidth',2)
plot3(lx(3:4),ly(3:4),lz(3:4),'red','LineWidth',2)
plot3(lx(5:6),ly(5:6),lz(5:6),'red','LineWidth',2)
plot3(lx(7:8),ly(7:8),lz(7:8),'red','LineWidth',2)
scatter3(lx(1:2),ly(1:2),lz(1:2),'filled','red')
scatter3(lx(3:4),ly(3:4),lz(3:4),'filled','red')
scatter3(lx(5:6),ly(5:6),lz(5:6),'filled','red')
scatter3(lx(7:8),ly(7:8),lz(7:8),'filled','red')
plot3([lx(1) lx(3)],[ly(1) ly(3)],[Emin Emin],'red','LineWidth',2)
plot3([lx(3) lx(5)],[ly(3) ly(5)],[Emin Emin],'red','LineWidth',2)
plot3([lx(5) lx(7)],[ly(5) ly(7)],[Emin Emin],'red','LineWidth',2)
plot3([lx(7) lx(1)],[ly(7) ly(1)],[Emin Emin],'red','LineWidth',2)
plot3([lx(1) lx(3)],[ly(1) ly(3)],[Emax Emax],'red','LineWidth',2)
plot3([lx(3) lx(5)],[ly(3) ly(5)],[Emax Emax],'red','LineWidth',2)
plot3([lx(5) lx(7)],[ly(5) ly(7)],[Emax Emax],'red','LineWidth',2)
plot3([lx(7) lx(1)],[ly(7) ly(1)],[Emax Emax],'red','LineWidth',2)
hold off
xlabel('k_x'), ylabel('k_y'), zlabel('E-E_F [eV]')
box on
grid on
axis([-2 2 -2 2 Emin Emax])
shading faceted
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

delta = 0.25;
for c=1:4
    t = text(kpts_(1,c),kpts_(2,c),Emax+delta,lbls{c});
    set(t,'FontSize',fs)
    set(t,'Color','red')
    set(t,'FontWeight','bold')
end
for c=1:4
    t = text(kpts_(1,c),kpts_(2,c),Emin-delta,lbls{c});
    set(t,'FontSize',fs)
    set(t,'Color','red')
    set(t,'FontWeight','bold')
end



%% fig1, plot trace bs
Efold = read_bin_mat([prefix 'Efold.bin']);
pos = read_bin_mat([prefix 'pos.bin']);
path = read_bin_mat([prefix 'path.bin']);

[Nb,Nk] = size(Efold);

% find indices for labels
J = zeros(1,5);
for c=2:4
    n = zeros(1,Nk);
    for d=1:Nk
        n(d) = norm(path(:,d)-kpts(:,c));
    end
    [~,j] = min(n);
    J(c) = j;
end
J(1) = 1; J(end) = Nk;

subplot(1,2,2), title('bandstructure trace')
hold on
for cb=I
    plot(pos,Efold(cb,:)-Ef,'LineWidth',2)
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


%% fig2, plot trace bs
figure

subplot(1,3,1), title('bandstructure')
hold on
for cb=I
    plot(pos,Efold(cb,:)-Ef,'LineWidth',2)
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


%% fig2, plot gradient bs
I = [vb_I(end-1:end)]; clr = {[0,0.447,0.741],[0.635,0.078,0.184]};
Gfold = read_bin_mat([prefix 'Gfold.bin']);
subplot(1,3,2), title('gradient, top 2 vb')
hold on
for cb=I
        
    g = squeeze(Gfold(cb,1,J(1):J(2)));
    scatter(pos(J(1):J(2)),g,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    g = squeeze(Gfold(cb,2,J(2):J(3)));
    scatter(pos(J(2):J(3)),g,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    g = squeeze(Gfold(cb,1,J(3):J(4)));
    scatter(pos(J(3):J(4)),g,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    g = squeeze(Gfold(cb,2,J(4):J(5)));
    scatter(pos(J(4):J(5)),g,'filled','MarkerFaceColor',clr{cb-I(1)+1});
end
hold off
set(gca,'XTick',pos(J));
set(gca,'XTickLabel',lbls);
box on
grid on
ylabel(['G [eV' char(197) ']'])
axis([pos(1) pos(end), -2.5 2.5])
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)


%% fig2, curvature bs
Cfold = read_bin_mat([prefix 'Cfold.bin']);
subplot(1,3,3), title('curvature, top 2 vb')
hold on
for cb=I
    c = squeeze(Cfold(cb,1,1,J(1):J(2)));
    scatter(pos(J(1):J(2)),c,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    c = squeeze(Cfold(cb,2,2,J(2):J(3)));
    scatter(pos(J(2):J(3)),c,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    c = squeeze(Cfold(cb,1,1,J(3):J(4)));
    scatter(pos(J(3):J(4)),c,'filled','MarkerFaceColor',clr{cb-I(1)+1});
    c = squeeze(Cfold(cb,2,2,J(4):J(5)));
    scatter(pos(J(4):J(5)),c,'filled','MarkerFaceColor',clr{cb-I(1)+1});
end
hold off
set(gca,'XTick',pos(J));
set(gca,'XTickLabel',lbls);
box on
grid on
ylabel(['C [eV' char(197) '^2]'])
axis([pos(1) pos(end), -15 15])
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

