% convert quaternions to euler angles in Bunge convention
function euler_angle=rod2euler(x)

[yaw, pitch, roll]=rod2angle(x,'ZXZ');
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
euler_angle=[yaw pitch roll];
