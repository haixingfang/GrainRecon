% General rotation: rotation angles of alpha, beta, gamma about z, y, x,
% respectively
function Rxyz=RotationMatrixFromAngle(gamma,beta,alpha)

Rz=[cosd(alpha) -sind(alpha) 0;sind(alpha) cosd(alpha) 0;0 0 1];
    
Ry=[cosd(beta) 0 sind(beta);0 1 0;-sind(beta) 0 cosd(beta)];

Rx=[1 0 0;0 cosd(gamma) -sind(gamma);0 sind(gamma) cosd(gamma)];

Rxyz=Rz*Ry*Rx; % yaw, pitch, raw



