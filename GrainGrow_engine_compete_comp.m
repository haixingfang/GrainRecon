% GrainGrow_engine_compete_comp_for_unindexed_voxels.m
% input indexed results
% for each un-indexed voxel, find OR candidates from a neiboring region
% calculate the corresponding completeness for each of the candidate OR
% take the one yielding the maximum completeness
% April 20, 2022
function [DS_out] = GrainGrow_engine_compete_comp(DS_in,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00, ...
                        P0y,P0z,RotAxisOffset,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,OutputFolder,search_radius)

if nargin<30
    search_radius=20; % radius of search region to find candidate orientations [pixel]
end
sprintf('search_radius = %d pixels',search_radius)

DS_in.SeedID = 1:max(DS_in.GrainId(:));
for i = DS_in.SeedID
    ind = ind2sub(size(DS_in.GrainId),find(DS_in.GrainId == i));
    if ~isempty(ind)
        DS_in.EulerAngle(i,:)=DS_in.EulerZXZ(:,ind(1))';
    end
end
DS_out = DS_in;

% indices for un-indexed voxels
un_indexed_voxels=[];
[un_indexed_voxels(:,1), un_indexed_voxels(:,2), un_indexed_voxels(:,3)]=ind2sub(size(DS_in.VisitFlag), ...
        find(DS_in.GrainId==0 & DS_in.Mask==1));
if ~isempty(un_indexed_voxels)
    sprintf('find %d un-indexed voxels among all %d voxels',length(un_indexed_voxels(:,1)),length(find(DS_in.Mask==1)))
else
    sprintf('find 0 un-indexed voxels')
end

nr_voxel_per_iter = 30000;
total_iter = floor(length(un_indexed_voxels(:,1))/nr_voxel_per_iter)+1;
sprintf('All un-indexed voxels divide to %d segments (each contains %d voxels) per gpu calculation.',total_iter,nr_voxel_per_iter)

for i=1:total_iter
    if i <= total_iter-1
        wait_revising_voxels = un_indexed_voxels(nr_voxel_per_iter*(i-1)+1:nr_voxel_per_iter*i,:);
    else
        wait_revising_voxels = un_indexed_voxels(nr_voxel_per_iter*(i-1)+1:end,:);
    end
    tic
    [id_neigb_all, id_neigb_ind] = find_neighbor_ORs(DS_in, wait_revising_voxels, search_radius);
    find_neighbor_time=toc;

    tic
    DS_out = grow_unindexed_compete_comp_cuda(DS_out,wait_revising_voxels,id_neigb_all,id_neigb_ind, ...
                        proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,RotDet,thetamax,lambda_min,lambda_max, ...
                        Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset,pixelysize,pixelzsize, ...
                        detysize,detzsize,BeamStopY,BeamStopZ,RecVolumePixel,tomo_scale,VoxSize,simap_data_flag);
    grow_time=toc;
    sprintf('The find_neighbor_time and the grow_time take %0.2f s and %0.2f s, respectively',find_neighbor_time,grow_time)
    sprintf('Iter %d: reconstructed volume fraction: %.4f = %d / %d ',i,length(find(DS_out.GrainId>0))/length(find(DS_out.Mask==1)), ...
        length(find(DS_out.GrainId>0)),length(find(DS_out.Mask==1)))

    save(fullfile(OutputFolder,strcat(num2str(i), 'DS.mat')),'DS_out','-v7.3');
end















