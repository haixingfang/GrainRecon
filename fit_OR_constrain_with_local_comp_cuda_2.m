% group all the sub gridded orientations around the index_seed OR
% forward calculate completeness by gpu_cuda
% select the orientation with the maximum completeness as initial values for fitting
% fitting once only
% 2 times of local gridding
% Updated on Oct 21, 2022

function [indexing_OR,confident_index_list,spots_pair]=fit_OR_constrain_with_local_comp_cuda_2(i,index_seeds,RotDet, ...
    proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
    thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ,fit_all_flag)

% By default, only select the first one for fitting
if nargin<31 || fit_all_flag~=1
    Nr_fit=1;
elseif fit_all_flag==1
%     Nr_fit=length(index_seeds(:,1));
    Nr_fit=10;
end

q_all=[];
for j=1:length(index_seeds(:,1))
%     OR_local=ori_local_sampling(index_seeds(j,5:8),10,0.001); % cover the whole range from -0.02 to 0.02, delta = nr * gridsize
    OR_local=ori_local_sampling(index_seeds(j,5:8),8,0.004);
    q_all((j-1)*length(OR_local.q(:,1))+1:j*length(OR_local.q(:,1)),1:4)=OR_local.q;
    q_all((j-1)*length(OR_local.q(:,1))+1:j*length(OR_local.q(:,1)),5)=j;
end

hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';

nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
pos_cuda = single(pos_indexing);
% q_cuda = single(reshape(q_all(:,1:4)',1,[]));
S_cuda = single(reshape(S',1,9));
B_cuda = single(reshape(B',1,9));
hkl_cuda = single(reshape(hkl,[1 size(hkl,1)*size(hkl,2)]));
proj_bin_bw_cuda = single(reshape(proj_bin_bw,1,[]));
rot_angles_cuda = single(rot_angles);
RotDet_cuda = single(reshape(RotDet',1,9));
param_cuda = single([nr_hkl nr_rot Lsam2sou Lsam2det P0y P0z dety00 detz00 ...
                    pixelysize pixelzsize detysize detzsize BeamStopY ...
                    BeamStopZ thetamax lambda_min lambda_max RotAxisOffset]);

if nr_rot * nr_hkl > 181*80
    error("error: the number of possible diffraction events exceeds the limit, \nplease increase the limit of the pre_allocated size for dis_simu_all in cuda_forward_comp.cu")
end

nr_ori_per_iter = 300000;
total_iter = floor(length(q_all(:,1))/nr_ori_per_iter)+1;
Output = zeros(size(q_all,1),9);
for ii=1:total_iter
    if ii <= total_iter-1
        wait_q = q_all(nr_ori_per_iter*(ii-1)+1:nr_ori_per_iter*ii,:);
    else
        wait_q = q_all(nr_ori_per_iter*(ii-1)+1:end,:);
    end
    q_cuda = single(reshape(wait_q(:,1:4)',1,[]));
    Output_cuda = cuda_forward_comp(pos_cuda, q_cuda, S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda);
    tmp = gather(Output_cuda);
    gpuDevice(1); % reset gpu device
    for j=1:size(wait_q,1)
        Output(j+(ii-1)*nr_ori_per_iter,1:4) = tmp((j-1)*4+1:(j-1)*4+4);
        Output(j+(ii-1)*nr_ori_per_iter,5:9) = wait_q(j,:);
    end
end

% do this multiple times
Output_max_temp=zeros(length(index_seeds(:,1)),9);
for j=1:length(index_seeds(:,1))
    if ~isempty(find(Output(:,9)==j))
        ind0=find(Output(:,3)==max(Output(Output(:,9)==j,3)) & Output(:,9)==j);
        Output_max_temp(j,:)=Output(ind0(1),:);
%         else
%             sprintf('Warning: OR %d returns empty result',j)
    end
end
Output_max_temp = sortrows(Output_max_temp,3,'descend');
calc_iter=3;
for ii=1:calc_iter
    [Output,Output_max_temp] = forward_calc_gridding_OR(Output_max_temp,5,0.001,pos_cuda,S_cuda, B_cuda, ...
                        hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda);
end

multi_fit=1; % fit more than 1 candidates?
if multi_fit==0
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
    [Nr_simu,Nr_intersect,dis_mean,~,SpotsPair]=forward_comp_pair_spots(UU,proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
    SpotsPair=double(SpotsPair);
    
    % save the result
    indexing_OR=[i pos_indexing Nr_simu Nr_intersect Nr_intersect/Nr_simu ...
        dis_mean yaw pitch roll rod quat];
    confident_index_list=Output(ind,9);
    spots_pair=[ones(size(SpotsPair,1),1)*i SpotsPair];
else
    Output_max=zeros(length(index_seeds(:,1)),9);
    for j=1:length(index_seeds(:,1))
        if ~isempty(find(Output(:,9)==j))
            ind0=find(Output(:,3)==max(Output(Output(:,9)==j,3)) & Output(:,9)==j);
            Output_max(j,:)=Output(ind0(1),:);
%         else
%             sprintf('Warning: OR %d returns empty result',j)
        end
    end
    Output_max = sortrows(Output_max,3,'descend');

    indexing_result=zeros(Nr_fit,18);
    for j=1:Nr_fit
        x0=double(Output_max(j,6:8));
        Err2Sum=@(x)fit_OR_fun_local(x,proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,RotDet, ...
            thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005}, ...
            'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10);
%         stop_fitting=0;
%         while(stop_fitting~=1)
            LB=[max([x0(1)-0.02 -1]) max([x0(2)-0.02 -1]) max([x0(3)-0.02 -1])];
            UB=[min([x0(1)+0.02 1]) min([x0(2)+0.02 1]) min([x0(3)+0.02 1])];
            [x,~,~,~] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
%             if sum(abs(x-LB)<0.001)>0 || sum(abs(x-UB)<0.001)>0
%                 x0=x;
%             else
%                 stop_fitting=1;
%             end
%         end
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
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        SpotsPair{j}=double(SpotsPair{j});
        indexing_result(j,:)=[i pos_indexing Nr_simu Nr_intersect Nr_intersect/Nr_simu ...
            dis_mean yaw pitch roll rod quat];
    end
    
    % save the result
    % the solution has the highest completeness value
    [~,ind]=max(indexing_result(:,7));
    indexing_OR=indexing_result(ind,:);
    confident_index_list=Output_max(ind,9);
    spots_pair=[ones(size(SpotsPair{ind},1),1)*i ...
                SpotsPair{ind}];
end
