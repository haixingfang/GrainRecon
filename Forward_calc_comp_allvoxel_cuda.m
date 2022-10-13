% call cuda_grow_indexed_region.cu
% April 12, 2022
function [Completeness_out,DisMedian_out,Ninter_out,Comp_max]=Forward_calc_comp_allvoxel_cuda(I,Imask,DisMedian,Ninter,voxel_waiting_growth0,maxD, ...
                    UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                    RotDet,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                    pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                    RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian)
% I : input image - completeness map

% % % for testing
% I = Completeness;
% Imask = Mask;
% DisMedian = Dismedian;

% %%%
Completeness_out = I;
DisMedian_out = DisMedian;
Ninter_out = Ninter;

voxel_waiting_growth = zeros(length(voxel_waiting_growth0(:,1)),6);
voxel_waiting_growth(:,1:3) = voxel_waiting_growth0;
for i=1:length(voxel_waiting_growth0(:,1))
    voxel_waiting_growth(i,4) = I(voxel_waiting_growth0(i,1),voxel_waiting_growth0(i,2),voxel_waiting_growth0(i,3));
    voxel_waiting_growth(i,5) = Imask(voxel_waiting_growth0(i,1),voxel_waiting_growth0(i,2),voxel_waiting_growth0(i,3));
    voxel_waiting_growth(i,6) = DisMedian(voxel_waiting_growth0(i,1),voxel_waiting_growth0(i,2),voxel_waiting_growth0(i,3));
end

hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';
nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
U_cuda = single(reshape(UU',1,9));
S_cuda = single(reshape(S',1,9));
B_cuda = single(reshape(B',1,9));

hkl_cuda = single(reshape(hkl,[1 size(hkl,1)*size(hkl,2)]));
proj_bin_bw_cuda = single(reshape(proj_bin_bw,1,[]));
rot_angles_cuda = single(rot_angles);
RotDet_cuda = single(reshape(RotDet',1,9));
param_cuda = single([nr_hkl nr_rot Lsam2sou Lsam2det P0y P0z ...
                    dety00 detz00 pixelysize pixelzsize detysize detzsize BeamStopY ...
                    BeamStopZ thetamax lambda_min lambda_max RotAxisOffset]);

% seed_comp=0;
% drop_off=0;
param1_cuda = single([RecVolumePixel(:,1)' simap_data_flag VoxSize tomo_scale.Center' tomo_scale.Dimension ...
                0 0 maxDmedian]);
if nr_rot * nr_hkl > 181*80
    error("error: the number of possible diffraction events exceeds the limit, \nplease increase the limit of the pre_allocated size for dis_simu_all in cuda_forward_comp.cu")
end

% change to 1D
voxel_waiting_growth = single(reshape(voxel_waiting_growth',1,[]));

% mexcuda cuda_grow_indexed_region.cu
clear Output_cuda tmp;
Output_cuda = cuda_grow_indexed_region(voxel_waiting_growth, U_cuda, ...
                        S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda, param1_cuda);
tmp = gather(Output_cuda);
gpuDevice(1);

% [xn yn zn Nr_simu Nr_intersect dis_median growth_flag]           
% note: the columes have different meanings from the cpu version
neg_list = zeros(length(voxel_waiting_growth)/6,7);
for i=1:length(voxel_waiting_growth)/6
    neg_list(i,1:7) = tmp((i-1)*7+1:(i-1)*7+7);
end
for j=1:length(neg_list(:,1))
    xn = neg_list(j,1);
    yn = neg_list(j,2);
    zn = neg_list(j,3);
    Completeness_out(xn,yn,zn) = neg_list(j,5)/neg_list(j,4);
    DisMedian_out(xn,yn,zn) = neg_list(j,6);
    Ninter_out(xn,yn,zn) = neg_list(j,5);
end
Comp_max=max(neg_list(:,5)./neg_list(:,4));

