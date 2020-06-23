% plots the primitive cell as well as the (rotated) super cell
% Mo: blue, S: green

prefix = 'output/l/';

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

%% primitive cell
close all
figure
hold on


[B,Ap,a,N,T,id,crystal] = read_poscar([prefix 'out.psc.ap']);
Ap = B*Ap; B(3,3) = 4.5;
I1 = T==1; I2 = T==2; clr = {'blue','green'};


pp = B*P1; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P2; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P3; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
scatter3(Ap(1,I1),Ap(2,I1),Ap(3,I1),'filled',clr{1},'SizeData',120)
scatter3(Ap(1,I2),Ap(2,I2),Ap(3,I2),'filled',clr{2},'SizeData',120)


[B,Ap,a,N,T,id,crystal] = read_poscar([prefix 'super.psc']);
Ap = B*Ap; B(:,3) = B(:,3)/norm(B(:,3))*4.5;
I1 = T==1; I2 = T==2; clr = {'blue','green'};

pp = B*P1; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P2; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
pp = B*P3; plot3(pp(1,:),pp(2,:),pp(3,:),'black','LineWidth',2)
scatter3(Ap(1,I1),Ap(2,I1),Ap(3,I1),'filled',clr{1},'SizeData',120)
scatter3(Ap(1,I2),Ap(2,I2),Ap(3,I2),'filled',clr{2},'SizeData',120)


hold off
axis equal
box on
grid on
xlabel('x'), ylabel('y'), zlabel('z')
set(gca,'FontSize',fs)
set(gca,'LineWidth',1.5)


