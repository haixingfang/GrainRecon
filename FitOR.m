% find a rotation matrix to transform euler_grains to euler_grains_ref
% x in radian
function Err2Sum=FitOR(x,euler_grains,euler_grains_ref,cs)

% x=x*180/pi; % [deg]
if true
% if false
    % using misori, sometimes it can be wrong when symmetry is not considered
    Su=euler2u(x(1)*pi/180,x(2)*pi/180,x(3)*pi/180);
    for j=1:length(euler_grains(:,1))
        u0=euler2u(euler_grains(j,1)*pi/180,euler_grains(j,2)*pi/180,euler_grains(j,3)*pi/180);
        u1=Su*u0;
        u_ref=euler2u(euler_grains_ref(j,1)*pi/180,euler_grains_ref(j,2)*pi/180,euler_grains_ref(j,3)*pi/180);
        [ang,~]=misori2(u1,u_ref);
%         Err(j)=sum(sum((u1-u_ref).^2));
%         Err(j)=sqrt(Err(j));
        Err(j)=ang;
    end
    Err2Sum=sum(Err)/length(euler_grains(:,1));
else
    % using mtex
    Su=euler2u(x(1)*pi/180,x(2)*pi/180,x(3)*pi/180);
    for j=1:length(euler_grains(:,1))
        u0=euler2u(euler_grains(j,1)*pi/180,euler_grains(j,2)*pi/180,euler_grains(j,3)*pi/180);
        u1=Su*u0;
        euler_grains_new=u2euler_corr(u1);
        rot_new = rotation('Euler',euler_grains_new(:,1)*degree,euler_grains_new(:,2)*degree,euler_grains_new(:,3)*degree);
        o_new = orientation(rot_new,cs);

        rot_ref = rotation('Euler',euler_grains_ref(:,1)*degree,euler_grains_ref(:,2)*degree,euler_grains_ref(:,3)*degree);
        o_ref = orientation(rot_ref,cs);
        ang=angle(o_new,o_ref(j))/degree;

        Err(j)=ang;    
    end
    Err2Sum=sum(Err)/length(euler_grains(:,1));
end
