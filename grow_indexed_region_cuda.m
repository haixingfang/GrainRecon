% call cuda_grow_indexed_region.cu
% April 12, 2022
function [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region_cuda(I,Imask,DisMedian,Ninter,pos_indices,Spot_L,drop_off,maxD, ...
        indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max, ...
        Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset, ...
        pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian)
        
%
% I : input image - completeness map
% J : logical output image of region
% pos_indices: [x,y,z] the position of the seedpoint
% drop_off : completeness decrease amount, by default 0.02
% minD : minimum acceptable completeness weighted center change [pixel]
%
% Given a bounding box whose dimensions are calculated according to the spot size for checking the voxels to grow or not
% March 30, 2022

% % for testing
% I = Completeness;
% Imask = Mask;
% DisMedian = Dismedian;
% i=1;
% Spot_L = mean(Spots_pair{i}(:,21))+3*std(Spots_pair{i}(:,21)); % bounding area of the spot [pixel^2]
% Spot_L = sqrt(Spot_L/pi); % bounding half length of the spot [pixel]

% %%%
x = pos_indices(1);
y = pos_indices(2);
z = pos_indices(3);
I(x,y,z) = indexed_comp;
J = zeros(size(I)); % Output
J(x,y,z) = 2;
Completeness_out = I;
DisMedian_out = DisMedian;
Ninter_out = Ninter;

seed_comp = I(x,y,z); % The seeding completeness of the segmented region
completeness_sum = seed_comp;
completeness_sum_pos = seed_comp*[x,y,z];
completeness_center0 = completeness_sum_pos/completeness_sum; % weighted center

Box_L = ceil(Spot_L*(mean([pixelysize pixelzsize])) / ((Lsam2det/Lsam2sou)*VoxSize)); % half width of the bounding box for the grain [voxel]
% define the boundary of the box
for i = 1:3
    if pos_indices(i)-Box_L > 0 
        BoxDim(i,1) = pos_indices(i)-Box_L;
    else
        BoxDim(i,1) = 1;
    end
    if pos_indices(i)+Box_L < size(I,i) 
        BoxDim(i,2) = pos_indices(i)+Box_L;
    else
        BoxDim(i,2) = size(I,i);
    end
end
[x_grid,y_grid,z_grid] = meshgrid(BoxDim(1,1):BoxDim(1,2),BoxDim(2,1):BoxDim(2,2),BoxDim(3,1):BoxDim(3,2));
voxel_waiting_growth0 = [reshape(x_grid,1,[])' reshape(y_grid,1,[])' reshape(z_grid,1,[])'];
% inner sphere instead of a cubic, added on Sep 12, 2022
CenterPoint = round(BoxDim(:,2)-BoxDim(:,1))/2+BoxDim(:,1);
CenterPoint = repmat(CenterPoint',[length(voxel_waiting_growth0(:,1)),1]);
Dis = sqrt((voxel_waiting_growth0(:,1) - CenterPoint(:,1)).^2 + ...
            (voxel_waiting_growth0(:,2) - CenterPoint(:,2)).^2 + ...
            (voxel_waiting_growth0(:,3) - CenterPoint(:,3)).^2);
voxel_waiting_growth0 = voxel_waiting_growth0(Dis<=Box_L,:);
% %%
sprintf('Investigate the bounding volume [%d %d; %d %d; %d %d] for the seeding voxel [%d %d %d] and check %d voxels for growth.', ...
    reshape(BoxDim',1,[]),pos_indices,length(voxel_waiting_growth0(:,1)))

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
param1_cuda = single([RecVolumePixel(:,1)' simap_data_flag VoxSize tomo_scale.Center' tomo_scale.Dimension ...
                seed_comp drop_off maxDmedian]);
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

% [xn yn zn Nr_simu Nr_intersect dis_median growth_flag]            
% note: the columes have different meanings from the cpu version
neg_list = zeros(length(voxel_waiting_growth)/6,7);
for i=1:length(voxel_waiting_growth)/6
    neg_list(i,1:7) = tmp((i-1)*7+1:(i-1)*7+7);
end

ind=find(neg_list(:,7)==1); % voxels to be accept
if ~isempty(ind)
    for j=1:length(ind)
        xn = neg_list(ind(j),1);
        yn = neg_list(ind(j),2);
        zn = neg_list(ind(j),3);
        J(xn,yn,zn)=2;
        Completeness_out(xn,yn,zn) = neg_list(ind(j),5)/neg_list(ind(j),4);
        DisMedian_out(xn,yn,zn) = neg_list(ind(j),6);
        Ninter_out(xn,yn,zn) = neg_list(ind(j),5);
    end
end

% Return the grown volume as logical matrix
J=J>1;
[indexed_indices(:,1),indexed_indices(:,2),indexed_indices(:,3)]=ind2sub(size(J),find(J>0));

if length(indexed_indices(:,1))>1
    for j=1:length(indexed_indices(:,1))
        completeness_sum = completeness_sum + ...
            Completeness_out(indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3));
        completeness_sum_pos = completeness_sum_pos + ...
            Completeness_out(indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3))* ...
            [indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3)];
    end
end

completeness_center = completeness_sum_pos/completeness_sum; % update weighted center
if sqrt(sum((completeness_center-completeness_center0).^2))>max([maxD maxD/10000*length(indexed_indices(:,1))]) % relative for pixel number > 10000
    center = round(completeness_center);
    replace_center = 1;
else
    center = completeness_center0;
    replace_center = 0;
end
gpuDevice(1);

