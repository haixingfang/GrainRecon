% June 7, 2021
% verified
function U = quaternion2U(x)

% x is represented by quaternions
% Extract the values from Q
q0 = x(1);
q1 = x(2);
q2 = x(3);
q3 = x(4);

% First row of the rotation matrix
r00 = 2 * (q0 * q0 + q1 * q1) - 1;
r01 = 2 * (q1 * q2 - q0 * q3);
r02 = 2 * (q1 * q3 + q0 * q2);

% Second row of the rotation matrix
r10 = 2 * (q1 * q2 + q0 * q3);
r11 = 2 * (q0 * q0 + q2 * q2) - 1;
r12 = 2 * (q2 * q3 - q0 * q1);

% Third row of the rotation matrix
r20 = 2 * (q1 * q3 - q0 * q2);
r21 = 2 * (q2 * q3 + q0 * q1);
r22 = 2 * (q0 * q0 + q3 * q3) - 1;

% 3x3 rotation matrix
U = [r00 r01 r02;r10 r11 r12;r20 r21 r22];

