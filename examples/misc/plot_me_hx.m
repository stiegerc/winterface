% figure 1
% Hamiltonian H(R) blocks such that only nn interactions exist.
% The corresponding R is given on top of each block
%
% figure 2
% Hamiltonian H(R) blocks where more than one neighbor in z exists.
% The corresponding R is given on top of each block
%
% figure 3
% Comparison of the bandstructures on the X-Gamma-X trace
% green: from nn interactions only
% blue: multiple R in z direction
% The expectation is:
% More green bands than blue ones, but all the blue ones have a green
% copy. This is due to the folding in z direction that is present for the
% green bands but absent for the blue ones. Deviations are due to the
% discarding of small interactions and the number of k-points used in the
% bandstructure calculations. Particularly for the nn interactions only
% version, bonds were filtered in x and z directions, whereas in the other
% case, filtering was done only along x. Accuracy can be improved by
% including more interactions at the cost of a dramatic increase in
% computation time.

prefix = 'output/hx/';
L = 1;   % depth of valence/conduction bands shown
fs = 20; % font size


close all

%% fig1
[H,R] = read_bin_sparse_hr('mos2_nn.hr');
NR = size(R,2);
N = ceil(sqrt(NR));

figure(1)
for c=1:NR
    subplot(N,N,c)
    spy(H{c})
    title([num2str(R(1,c)) ' ' num2str(R(2,c)) ' ' num2str(R(3,c))])
end

%% fig2
[H,R] = read_bin_sparse_hr('mos2.hr');
NR = size(R,2);
N = ceil(sqrt(NR));

figure(2)
for c=1:NR
    subplot(N,N,c)
    spy(H{c})
    title([num2str(R(1,c)) ' ' num2str(R(2,c)) ' ' num2str(R(3,c))])
end

%% fig3
Enn = read_bin_mat([prefix 'nn_Escal.bin']);
E = read_bin_mat([prefix 'Escal.bin']);
Ef = read_outcar('data/OUTCAR');

vb = max(Enn(Enn<Ef));
cb = min(Enn(Enn>Ef));

vbinn = find(Enn(:,1)<Ef,1,'last');
vbi = find(E(:,1)<Ef,1,'last');

Nk = size(E,2);
kpl = linspace(-.5,.5,Nk);

figure(3)

subplot(1,2,1)
hold on
plot(kpl,Enn(1:vbinn,:)-Ef,'green','LineWidth',3)
plot(kpl,E(1:vbi,:)-Ef,'blue','LineWidth',1.5)
hold off
axis([-.5 .5 vb-Ef-L vb-Ef+.1])
title('valence bands')
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel('E-E_F [eV]')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)


subplot(1,2,2)
hold on
plot(kpl,Enn(vbinn+1:end,:)-Ef,'green','LineWidth',3)
plot(kpl,E(vbi+1:end,:)-Ef,'blue','LineWidth',1.5)
hold off
axis([-.5 .5 cb-Ef-.1 cb-Ef+L])
title('conduction bands')
set(gca,'XTick',[-.5 0 .5]);
set(gca,'XTickLabel',{'X','\Gamma','X'});
box on
grid on
ylabel('E-E_F [eV]')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

