function RotDet=get_det_R(tilt_x,tilt_y,tilt_z)

fprintf('tilt angles = [%.4f %.4f %.4f] (deg)\n',tilt_x,tilt_y,tilt_z);
RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
RotDet=RotZ*RotY*RotX;