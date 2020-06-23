interf = wbh('wbhs/interf.wbh');
dbl_left = wbh('wbhs/dbl_left.wbh');
dbl_right = wbh('wbhs/dbl_right.wbh');



%% edit wads

% rotate to x axis
phi = -90; phi = phi/180*pi;
R = [
    cos(phi) -sin(phi) 0
    sin(phi) cos(phi) 0
    0 0 1
    ];

% swap periodic axis z,y
phi = 90; phi = phi/180*pi;
Rx = [
    1 0 0
    0 cos(phi) -sin(phi)
    0 sin(phi) cos(phi)
    ];

R = Rx*R;

interf.cell.B = R*interf.cell.B;
interf.write('wbhs/interf_rot.wbh');

dbl_left.cell.B = R*dbl_left.cell.B;
dbl_left.cell.B(1,:) = -dbl_left.cell.B(1,:); % invert x direction
dbl_left.write('wbhs/dbl_left_rot.wbh');

dbl_right.cell.B = R*dbl_right.cell.B;
dbl_right.cell.B(:,3) = -dbl_right.cell.B(:,3); % invert vacuum direction
dbl_right.write('wbhs/dbl_right_rot.wbh')


