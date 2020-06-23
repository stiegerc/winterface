% super cell in black, primitive cell in red
% Mo: blue, S: green

prefix = 'output/p/';
[sB,sAp,sa,sN,sT,sid,scrystal] = read_poscar([prefix 'super.psc']);
sAp = sB*sAp;
n = [norm(sB(:,1)) norm(sB(:,2)) norm(sB(:,3))];
[~,j] = max(n); sB(:,j) = 4.5*sB(:,j)/norm(sB(:,j));

[pB,pAp,pa,pN,pT,pid,pcrystal] = read_poscar([prefix 'prim.psc']);
pAp = pB*pAp;
n = [norm(pB(:,1)) norm(pB(:,2)) norm(pB(:,3))];
[~,j] = max(n); pB(:,j) = 4.5*pB(:,j)/norm(pB(:,j));

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


close all
figure

hold on
pp = sB*P1; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = sB*P2; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = sB*P3; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
I1 = sT==1; I2 = sT==2; clr = {'blue','green'};
scatter3(sAp(1,I1),sAp(2,I1),sAp(3,I1),'filled',clr{1},'SizeData',120)
scatter3(sAp(1,I2),sAp(2,I2),sAp(3,I2),'filled',clr{2},'SizeData',120)

pp = pB*P1; plot3(pp(1,:),pp(2,:),pp(3,:),'red','LineWidth',2)
pp = pB*P2; plot3(pp(1,:),pp(2,:),pp(3,:),'red','LineWidth',2)
pp = pB*P3; plot3(pp(1,:),pp(2,:),pp(3,:),'red','LineWidth',2)
I1 = pT==1; I2 = pT==2; clr = {'blue','green'};
scatter3(pAp(1,I1),pAp(2,I1),pAp(3,I1),'filled',clr{1},'SizeData',120)
scatter3(pAp(1,I2),pAp(2,I2),pAp(3,I2),'filled',clr{2},'SizeData',120)
hold off

axis equal
box on
grid on
xlabel('x'), ylabel('y'), zlabel('z')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)

