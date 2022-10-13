% first forward calculations for locally sampled orientations 
% them, fit to obtain the precise orientation
% September 27, 2021
% speed test: 12.4 times slower than fit_OR_constrain_with_local_comp.m

function [indexing_OR,indexing_result,confident_index_list,spots_pair]=fit_OR_constrain_with_local_comp_mex(i,index_seeds,RotDet, ...
    proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
    thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
    RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

SpotsPair={};
indexing_result=zeros(length(index_seeds(:,1)),18);
if length(index_seeds(:,1))>20
    parfor j=1:length(index_seeds(:,1))
        OR_local=ori_local_sampling(index_seeds(j,5:8),5);
        Output = ForwardCalc_OR_space_comp_mex(OR_local.q_b,OR_local.q_c,OR_local.q_d,size(OR_local.q_b,1),RotDet, ...
                    proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                    RotAxisOffset,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ);
        [~,ind]=max(Output(:,3));
        x0=double(Output(ind,6:8));
        Err2Sum=@(x)fit_OR_fun_local(x,proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,RotDet, ...
            thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005}, ...
            'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10);
        LB=[max([x0(1)-0.02 -1]) max([x0(2)-0.02 -1]) max([x0(3)-0.02 -1])];
        UB=[min([x0(1)+0.02 1]) min([x0(2)+0.02 1]) min([x0(3)+0.02 1])];
        [x,~,~,~] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
    
        [yaw, pitch, roll]=rod2angle(quat2rod([sqrt(1-x(1)^2-x(2)^2-x(3)^2) x(1) x(2) x(3)]),'ZXZ');
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
        rod=angle2rod(yaw*pi/180,pitch*pi/180,roll*pi/180,'ZXZ'); % Rodrigues
        quat=rod2quat(rod); % quaternion
        UU=quaternion2U(quat);
    
        % indexing verification
        [Nr_simu,Nr_intersect,dis_mean,~,SpotsPair{j}]=forward_comp_pair_spots(UU,proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        SpotsPair{j}=double(SpotsPair{j});
        % save the result
        indexing_result(j,:)=[i pos_indexing Nr_simu Nr_intersect Nr_intersect/Nr_simu ...
            dis_mean yaw pitch roll rod quat];
    end
else
    for j=1:length(index_seeds(:,1))      
        OR_local=ori_local_sampling(index_seeds(j,5:8),5);
        Output = ForwardCalc_OR_space_comp_mex(OR_local.q_b,OR_local.q_c,OR_local.q_d,size(OR_local.q_b,1),RotDet, ...
                    proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                    pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ);
        [~,ind]=max(Output(:,3));
        x0=double(Output(ind,6:8));
        Err2Sum=@(x)fit_OR_fun_local(x,proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,RotDet, ...
            thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005}, ...
            'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10);
        LB=[max([x0(1)-0.02 -1]) max([x0(2)-0.02 -1]) max([x0(3)-0.02 -1])];
        UB=[min([x0(1)+0.02 1]) min([x0(2)+0.02 1]) min([x0(3)+0.02 1])];
        [x,~,~,~] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
    
        [yaw, pitch, roll]=rod2angle(quat2rod([sqrt(1-x(1)^2-x(2)^2-x(3)^2) x(1) x(2) x(3)]),'ZXZ');
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
        rod=angle2rod(yaw*pi/180,pitch*pi/180,roll*pi/180,'ZXZ'); % Rodrigues
        quat=rod2quat(rod); % quaternion
        UU=quaternion2U(quat);
    
        % indexing verification
        [Nr_simu,Nr_intersect,dis_mean,~,SpotsPair{j}]=forward_comp_pair_spots(UU,proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        SpotsPair{j}=double(SpotsPair{j});
        % save the result
        indexing_result(j,:)=[i pos_indexing Nr_simu Nr_intersect Nr_intersect/Nr_simu ...
            dis_mean yaw pitch roll rod quat];
    end
end
% the solution has the highest completeness value
[~,confident_index_list]=max(indexing_result(:,7));
indexing_OR=indexing_result(confident_index_list,:);
spots_pair=[ones(size(SpotsPair{confident_index_list},1),1)*i ...
    SpotsPair{confident_index_list}];
end


function Output = ForwardCalc_OR_space_comp_mex(q_b,q_c,q_d,dim_OR,RotDet, ...
                    proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                    pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ)

% % for testing
% q_b=OR.q_b;
% q_c=OR.q_c;
% q_d=OR.q_d;
% dim_OR=size(OR.q_b,1);
% pos=pos_indexing;

    L=Lsam2sou+Lsam2det; % [mm]
    parfor i=1:dim_OR
        x = [q_b(i) q_c(i) q_d(i)];
        % x is represented by the last three quantities of the quaternions (a, b, c, d)
        % Extract the values from Q
        if (x(1)^2+x(2)^2+x(3)^2)<=1
            q0 = sqrt(1-(x(1)^2+x(2)^2+x(3)^2));
            q1 = x(1);
            q2 = x(2);
            q3 = x(3);

            % First row of the rotation matrix
            r00 = 2 * (q0 * q0 + q1 * q1) - 1;
            r01 = 2 * (q1 * q2 - q0 * q3);
            r02 = 2 * (q1 * q3 + q0 * q2);

            % Second row of the rotation matrix
            r10 = 2 * (q1 * q2 + q0 * q3);
            r11 = 2 * (q0 * q0 + q2 * q2) - 1;
            r12 = 2 * (q2 * q3 - q0 * q1);

            % Third row of the rotation matrix
            r20 = 2 * (q1 * q3 - q0 * q2);
            r21 = 2 * (q2 * q3 + q0 * q1);
            r22 = 2 * (q0 * q0 + q3 * q3) - 1;

            % 3x3 rotation matrix
            U = [r00 r01 r02;r10 r11 r12;r20 r21 r22];
            [Nr_simu,Nr_intersect,dis_median] = forward_comp_mex(pos,U,S,B,Ahkl,nrhkl, ...
                                    proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                                    P0y,P0z,RotAxisOffset,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                                    BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max);
            Output(i,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median q0 q1 q2 q3];
        end
        i;
    end
end
