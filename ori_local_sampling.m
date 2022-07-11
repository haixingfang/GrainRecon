% sampling around a specified orientation in a local orientation space
% Aug 24, 2021
% ori_center is expressed as quaternion or
% create (nr+1)^3 orientations

function OR_local=ori_local_sampling(ori_center,nr,gridsize)

% % for testing
% ori_center=quat;
% nr=5;

if length(ori_center)==3
    ori0(2:4)=ori_center;
    ori0(1)=sqrt(1-(ori0(2)^2+ori0(3)^2+ori0(4)^2));
else
    ori0=ori_center;
end
if nargin<3
    gridsize=0.004; % 0.0005 for 0.06 deg;0.001 for 0.12 deg; 0.005 for 0.60 deg
end

bgrid=ori0(2)-gridsize*nr/2:gridsize:ori0(2)+gridsize*nr/2;
cgrid=ori0(3)-gridsize*nr/2:gridsize:ori0(3)+gridsize*nr/2;
dgrid=ori0(4)-gridsize*nr/2:gridsize:ori0(4)+gridsize*nr/2;

q=ori0;
for i=1:length(bgrid)
    for j=1:length(cgrid)
        for k=1:length(dgrid)
            q=[q;sqrt(1-(bgrid(i)^2+cgrid(j)^2+dgrid(k)^2)) bgrid(i) cgrid(j) dgrid(k)];
        end
    end
end
q=unique(q,'rows');
rod=quat2rod(q);
% q=rod2quat(rod);
q_b=q(:,2);
q_c=q(:,3);
q_d=q(:,4);
[euler_angles(:,1), euler_angles(:,2), euler_angles(:,3)] = rod2angle(rod,'ZXZ');
euler_angles=euler_angles*180/pi;
for i=1:length(euler_angles(:,1))
    if euler_angles(i,1)<0
        euler_angles(i,1)=euler_angles(i,1)+360;
    end
    if euler_angles(i,2)<0
        euler_angles(i,2)=euler_angles(i,2)+360;
    end
    if euler_angles(i,3)<0
        euler_angles(i,3)=euler_angles(i,3)+360;
    end
end
OR_local.len=length(q_b(:,1));
OR_local.euler_angles=euler_angles;
OR_local.rod=rod;
OR_local.q=q;
OR_local.q_b=q_b;
OR_local.q_c=q_c;
OR_local.q_d=q_d;
    
% randnumber=randi([0 27000],1);
% quat=OR.q(randnumber,:);
% quat1=quat;
% quat2=quat1;
% quat2(4)=quat1(4)-0.001;
% % quat2(2:4)=quat1(2:4)-0.005;
% quat2(1)=sqrt(1-(quat2(2)^2+quat2(3)^2+quat2(4)^2));
% U1=quaternion2U(quat1);
% U2=quaternion2U(quat2);
% misori2(U1,U2)




