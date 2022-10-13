function [RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru)
% det_normal: along the beam
% detdiru: pointing to the right of the detector (opposite lab-y)

tilt_x = acosd(dot(detdiru,[0 0 1]))-90;
tilt_y = acosd(dot(det_normal.*[1 0 1],[0 0 1]))-90;
tilt_z = 90-acosd(dot(det_normal.*[1 1 0],[0 1 0]));
% fprintf('tilt angles = [%.4f %.4f %.4f] (deg)\n',tilt_x,tilt_y,tilt_z);
RotDet=get_det_R(tilt_x,tilt_y,tilt_z);
