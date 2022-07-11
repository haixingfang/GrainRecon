% divide orientation space using mtex
function OR=ori_division_mtex(crystalSymmetry_name,angular_resolution)

if nargin<2
    angular_resolution='2.7';
end
% define symmetries
% cs = crystalSymmetry('cubic');
cs = crystalSymmetry(crystalSymmetry_name);
grid_flag=0; % 0-based on mtex; 1-orthogonal gridding
if grid_flag==1
    gridsize=0.01;
    Q_b=-0.38:gridsize:0.38;
    Q_c=-0.38:gridsize:0.38;
    Q_d=-0.38:gridsize:0.35;
    count=0;
    for i=1:length(Q_b)
        for j=1:length(Q_c)
            for k=1:length(Q_d)
                count=count+1;
                q_b(count,1)=Q_b(i);
                q_c(count,1)=Q_c(j);
                q_d(count,1)=Q_d(k);
            end
        end
    end
    q_a=sqrt(1-(q_b.^2+q_c.^2+q_d.^2));
%         [euler_angles1,euler_angles2,euler_angles3]=rod2angle(quat2rod([q_a q_b q_c q_d]),'ZXZ');
    r = rotation(quaternion(q_a,q_b,q_c,q_d));
    ori=orientation(r,cs)
else
    % define a grid of orientations
%     ori = equispacedSO3Grid(cs,'resolution',1*degree)
    if strcmp(angular_resolution,'0.5')
        ori = equispacedSO3Grid(cs,'resolution',0.5*degree)
    elseif strcmp(angular_resolution,'1')
        ori = equispacedSO3Grid(cs,'resolution',1*degree)
    elseif strcmp(angular_resolution,'2')
        ori = equispacedSO3Grid(cs,'resolution',2*degree)
    else
        ori = equispacedSO3Grid(cs,'resolution',2.7*degree)
    end
    if size(ori.phi1,1)>1
        euler_angles(:,1)=reshape(ori.phi1,size(ori.phi1,1)*size(ori.phi1,2),1);
        euler_angles(:,2)=reshape(ori.Phi,size(ori.Phi,1)*size(ori.Phi,2),1);
        euler_angles(:,3)=reshape(ori.phi2,size(ori.phi2,1)*size(ori.phi2,2),1);  
    else
        euler_angles(:,1)=ori.phi1';
        euler_angles(:,2)=ori.Phi';
        euler_angles(:,3)=ori.phi2';
    end
    euler_angles=euler_angles*180/pi;
    if size(ori.phi1,1)==1
        rod=angle2rod(ori.phi1',ori.Phi',ori.phi2','ZXZ');
    else
        ori_phi1=reshape(ori.phi1,1,[]);
        ori_Phi=reshape(ori.Phi,1,[]);
        ori_phi2=reshape(ori.phi2,1,[]);
        rod=angle2rod(ori_phi1',ori_Phi',ori_phi2','ZXZ');
    end
    q=rod2quat(rod);
    q_b=q(:,2);
    q_c=q(:,3);
    q_d=q(:,4);
end
ori_nr=length(q_b(:,1));

OR.len=length(q_b(:,1));
OR.euler_angles=euler_angles;
OR.rod=rod;
OR.q=q;
OR.q_b=q_b;
OR.q_c=q_c;
OR.q_d=q_d;







