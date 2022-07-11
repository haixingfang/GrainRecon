% call cuda_forward_comp_compete.cu
% April 20, 2022
function [DS_out] = grow_unindexed_compete_comp_cuda(DS_in,wait_revising_voxels,id_neigb_all,id_neigb_ind, ...
                    proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                    RotDet,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                    RecVolumePixel,tomo_scale,VoxSize,simap_data_flag)

DS_out = DS_in;

voxel_waiting_growth = [wait_revising_voxels id_neigb_ind]; % merge

hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';
nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
id_neigb_all_cuda = single(id_neigb_all);
Euler_cuda = single(reshape(DS_in.EulerAngle',1,[]));
S_cuda = single(reshape(S',1,9));
B_cuda = single(reshape(B',1,9));

hkl_cuda = single(reshape(hkl,[1 size(hkl,1)*size(hkl,2)]));
proj_bin_bw_cuda = single(reshape(proj_bin_bw,1,[]));
rot_angles_cuda = single(rot_angles);
RotDet_cuda = single(reshape(RotDet',1,9));
param_cuda = single([nr_hkl nr_rot Lsam2sou Lsam2det P0y P0z ...
                    dety00 detz00 pixelysize pixelzsize detysize detzsize BeamStopY ...
                    BeamStopZ thetamax lambda_min lambda_max]);
param1_cuda = single([RecVolumePixel(:,1)' simap_data_flag VoxSize tomo_scale.Center' tomo_scale.Dimension]);
if nr_rot * nr_hkl > 181*80
    error("error: the number of possible diffraction events exceeds the limit, \nplease increase the limit of the pre_allocated size for dis_simu_all in cuda_forward_comp.cu")
end

% change to 1D
voxel_waiting_growth = single(reshape(voxel_waiting_growth',1,[]));

% mexcuda cuda_forward_comp_compete.cu
clear Output_cuda tmp;
Output_cuda = cuda_forward_comp_compete(voxel_waiting_growth, id_neigb_all_cuda, ...
                        Euler_cuda, S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda, param1_cuda);
tmp = gather(Output_cuda);

% [idx phi1 PHI phi2 GrainId Completeness Nr_intersect dis_median]
grow_list = zeros(length(voxel_waiting_growth)/5,8);
for i=1:length(voxel_waiting_growth)/5
    grow_list(i,:) = tmp((i-1)*8+1:(i-1)*8+8);
end

% update results
for i=1:length(wait_revising_voxels(:,1))
    xn = wait_revising_voxels(i,1);
    yn = wait_revising_voxels(i,2);
    zn = wait_revising_voxels(i,3);
    DS_out.Completeness(xn,yn,zn) = grow_list(i,6);
    DS_out.GrainId(xn,yn,zn) = grow_list(i,5);
    DS_out.Dismedian(xn,yn,zn) = grow_list(i,8);
    DS_out.Ninter(xn,yn,zn) = grow_list(i,7);
    DS_out.VisitFlag(xn,yn,zn) = 0.5;

    DS_out.EulerZXZ(:,xn,yn,zn) = grow_list(i,2:4)';
    DS_out.Rodrigues(:,xn,yn,zn) = angle2rod(grow_list(i,2)*pi/180,grow_list(i,3)*pi/180,grow_list(i,4)*pi/180,'ZXZ')';
    DS_out.IPF001(:,xn,yn,zn) = grow_list(i,2:4)'./norm(grow_list(i,2:4));
end
gpuDevice(1);

%{
% verification
i=2;
xn=wait_revising_voxels(i,1);
yn=wait_revising_voxels(i,2);
zn=wait_revising_voxels(i,3);
pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*tomo_scale.VoxSize(1)+tomo_scale.Center'; % [mm]
if simap_data_flag==1
    pos(1)=-pos(1);
    pos(2)=-pos(2);
end
id_neigb = id_neigb_all(id_neigb_ind(i,1):id_neigb_ind(i,2)); % neighbors ID

Comp_all=[];
dis_median_all=[];
Nr_intersect_all=[];
for j=1:length(id_neigb)
    UU=euler2u(DS_in.EulerAngle(id_neigb(j),1)*pi/180,DS_in.EulerAngle(id_neigb(j),2)*pi/180,DS_in.EulerAngle(id_neigb(j),3)*pi/180);
    [Nr_simu,Nr_intersect,dis_median]=forward_comp(UU,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
    Comp_all=[Comp_all;Nr_intersect/Nr_simu];
    dis_median_all=[dis_median_all;dis_median];
    Nr_intersect_all=[Nr_intersect_all;Nr_intersect];
end
[~, ind] = max(Comp_all);
Output = [i DS_in.EulerAngle(id_neigb(ind),:) id_neigb(ind) Comp_all(ind) Nr_intersect_all(ind) dis_median_all(ind)]
%}
