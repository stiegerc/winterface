% top row: bandstructure along the X-Gamma-X trace
% bottom row: first derivative of the top row in direction of the trace

fs = 20; % font size
Nvb = 10; % number of valence bands to plot
Ncb = 10; % number of conduction bands to plot

prefix = 'output/bl/';

[~,~,~,~,~,~,~,~,~,cb,vb] = ...
    read_omen_input('../mos2/ph_mat_par','../mos2/lattice_dat');
Ef = (vb+cb)/2; % fermi energy set to the middle of the gap

E = read_bin_mat([prefix 'Escal.bin']);
G = read_bin_mat([prefix 'Gscal.bin']);
C = read_bin_mat([prefix 'Cscal.bin']);

[Nb,Nk] = size(E);
e = E(:,1); vbi = find(e<Ef,1,'last');

vbI = 1:vbi;
cbI = vbi+1:Nb;

close all
figure
kpl = linspace(-.5,.5,Nk);

subplot(2,2,1)
hold on
plot(kpl,E(cbI(1:Ncb),:)-Ef,'LineWidth',3)
hold off
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel('E-E_F [eV]')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)
title('conduction bands')

subplot(2,2,2)
plot(kpl,E(vbI(end-Nvb+1:end),:)-Ef,'LineWidth',3)
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel('E-E_F [eV]')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)
title('valence bands')

ax = gca;
CC = ax.ColorOrder;
Nc = size(CC,1);


subplot(2,2,3)
hold on
g = squeeze(G(cbI(1:Ncb),1,:));
for c=1:size(g,1)
    scatter(kpl,g(c,:),'filled','SizeData',15,'MarkerFaceColor',CC(mod(c-1,Nc)+1,:))
end
hold off
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel(['G [eV' char(197) ']'])
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)


subplot(2,2,4)
hold on
g = squeeze(G(vbI(end-Nvb+1:end),1,:));
for c=1:size(g,1)
    scatter(kpl,g(c,:),'filled','SizeData',15,'MarkerFaceColor',CC(mod(c-1,Nc)+1,:))
end
hold off
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel(['G [eV' char(197) ']'])
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

