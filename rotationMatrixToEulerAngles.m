function [xyz]=rotationMatrixToEulerAngles(R)
Rt=R';
if det(eye(3)-R*Rt)<1e-6
    sy=sqrt(R(1,1)*R(1,1)+R(2,1)*R(2,1));
    if sy<1e-6
        x = atan2(-R(2,3),R(2,2));
        y = atan2(-R(3,1),sy);
        z = 0;
    else
        x = atan2(R(3,2),R(3,3));
        y = atan2(-R(3,1),sy);
        z = atan2(R(2,1),R(1,1));
    end
    xyz=[x y z]*180/pi;
else
    xyz=[];
end