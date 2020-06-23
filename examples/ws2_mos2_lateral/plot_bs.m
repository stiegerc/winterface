% This will plot the band edges of the local bandstructures (with kx = [-.5,...,.5]),
% where the cells were sorted along the transport direction.
% The switch in material is visible in the middle and on each side.
% The borders of the DFT simulation domain are given by dashed horizontal
% lines
% plotting the bandstructures themselves is optional as it is hard to see
% anything. 1:10 are from WS_2, 11:20 from MoS_2

plot_all_bs = true;
I = 1:20;

%% read data
Ap = cell(1,numel(I));
id = cell(1,numel(I));
E = cell(1,numel(I));

for j=I
    [~,Ap{j},~,~,~,id{j}] = read_poscar(['bsout/local_' sprintf('%03d',j) '.psc']);
    E{j} = read_bin_mat(['bsout/E_trace_local_' sprintf('%03d',j) '.bin']);
end
k = read_bin_mat('bsout/path.bin');
Nk = size(k,2); Nb = size(E{1},1);

Ef = read_outcar('data/OUTCAR');



%% atoms in order along transport is, sorting manually
% 'W-001','W-006','W-002','W-007','W-003','W-008','W-004','W-009','W-005','W-010'
% 'Mo-001','Mo-006','Mo-002','Mo-007','Mo-003','Mo-008','Mo-004','Mo-009','Mo-005','Mo-010'
E_ = cell(1,numel(I));
E_{20} = E{1}; E_{18} = E{2}; E_{16} = E{3}; E_{14} = E{4}; E_{12} = E{5};
E_{10} = E{6}; E_{8} = E{7}; E_{6} = E{8}; E_{4} = E{9}; E_{2} = E{10};
E_{19} = E{11}; E_{17} = E{12}; E_{15} = E{13}; E_{13} = E{14}; E_{11} = E{15};
E_{9} = E{16}; E_{7} = E{17}; E_{5} = E{18}; E_{3} = E{19}; E_{1} = E{20};
E = E_;


%% get band edges
vb = zeros(1,numel(I));
cb = zeros(1,numel(I));
for j=I
    vb(j) = max(E{j}(E{j}(:)<Ef));
    cb(j) = min(E{j}(E{j}(:)>Ef));
end

%% wrap around on each side
x = [-2 -1 0 I max(I)+1 max(I)+2 max(I)+3];
vb = [vb(end-2:end) vb vb(1:3)];
cb = [cb(end-2:end) cb cb(1:3)];

%% plot band egdes
close all

figure(1)
hold on
plot(x,vb-Ef,'blue'), scatter(x,vb-Ef,'filled','blue')
plot(x,cb-Ef,'green'), scatter(x,cb-Ef,'filled','green')
plot([1 1],[-.8 2],'k--')
plot([20 20],[-.8 2],'k--')
hold off
axis([min(x) max(x) -.8 2])
xlabel('cell index')
ylabel('E-E_F [eV]')



%% plot bandstructures
if plot_all_bs
        figure(2)

        kpl = ones(Nb,1)*k(1,:);
        for j=I
            h = subplot(4,5,j);
            scatter(h,kpl(:),E{j}(:),'.','black')
            title(h,num2str(j))
        end
end

