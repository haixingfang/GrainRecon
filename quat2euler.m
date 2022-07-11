% convert quaternions to euler angles in Bunge convention
function [yaw, pitch, roll]=quat2euler(x)

if length(x)==3
    [yaw, pitch, roll]=rod2angle(quat2rod([sqrt(1-x(1)^2-x(2)^2-x(3)^2) x(1) x(2) x(3)]),'ZXZ');
else
    [yaw, pitch, roll]=rod2angle(quat2rod([x(1) x(2) x(3) x(4)]),'ZXZ');
end

yaw=yaw*180/pi;
if yaw<0
    yaw=yaw+360;
end
pitch=pitch*180/pi;
if pitch<0
    pitch=pitch+360;
end
roll=roll*180/pi;
if roll<0
    roll=roll+360;
end
