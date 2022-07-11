function [q,rod]=ang2quat(euler_angles,unit_name)
if nargin>1
    if unit_name~='rad'
        euler_angles=euler_angles*pi/180;
    end
else
    euler_angles=euler_angles*pi/180;
end
rod=angle2rod(euler_angles(:,1),euler_angles(:,2),euler_angles(:,3),'ZXZ');
q=rod2quat(rod);