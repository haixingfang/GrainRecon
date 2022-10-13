% GrainGrow_engine_compete_comp_for_unindexed_voxels.m
% input indexed results
% for each un-indexed voxel, find OR candidates from a neiboring region
% calculate the corresponding completeness for each of the candidate OR
% take the one yielding the maximum completeness
% April 20, 2022
function [DS_out] = revise_all_unindexed_voxel_cuda(DS_in,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00, ...
                        P0y,P0z,RotAxisOffset,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag, ...
                        correct_already_indexed_only,search_radius)

if nargin<30 || ~exist('search_radius','var')
    search_radius=20; % radius of search region to find candidate orientations [pixel]
end
sprintf('search_radius = %d pixels',search_radius)

% DS_in.SeedID = 1:max(DS_in.GrainId(:));
% for i = DS_in.SeedID
%     ind = ind2sub(size(DS_in.GrainId),find(DS_in.GrainId == i));
%     if ~isempty(ind)
%         DS_in.EulerAngle(i,:)=DS_in.EulerZXZ(:,ind(1))';
%     end
% end
DS_in.SeedID = 1:max(DS_in.GrainId(:));
for i = DS_in.SeedID
    [x,y,z] = ind2sub(size(DS_in.GrainId),find(DS_in.GrainId == i));
    X = mean(x);
    Y = mean(y);
    Z = mean(z);
    DS_in.Coord(i,:) = [X,Y,Z]; % each row represents coodinate for each grain
    DS_in.nVox(i,1) = length(find(DS_in.GrainId == i)); % number of voxels for each grain
    DS_in.EulerAngle(i,:)=DS_in.EulerZXZ(:,x(1),y(1),z(1))';

    pos=((DS_in.Coord(i,:)+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*tomo_scale.VoxSize'+tomo_scale.Center'; % [mm]
    if simap_data_flag==1
        pos(1:2)=-pos(1:2);
    end
    U=euler2u(DS_in.EulerAngle(i,1)*pi/180,DS_in.EulerAngle(i,2)*pi/180,DS_in.EulerAngle(i,3)*pi/180);
    [Nr_simu,Nr_intersect,~]=forward_comp(U,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
    DS_in.SeedComp(i,1)=Nr_intersect/Nr_simu;
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

% add on June 7, 2022, correct too low completeness data
correct_low_C=1;
doubt_indexed_voxels=[];
[doubt_indexed_voxels(:,1),doubt_indexed_voxels(:,2),doubt_indexed_voxels(:,3)]=ind2sub(size(DS_in.GrainId), ...
        find(DS_in.Completeness>0 & DS_in.Completeness<0.4));
if ~isempty(doubt_indexed_voxels)
    sprintf('find %d low-C (smaller than 0.40) indexed voxels among all %d voxels',length(doubt_indexed_voxels(:,1)),length(find(DS_in.Mask==1)))
else
    sprintf('find 0 low-C indexed voxels')
end
if correct_low_C==1
    un_indexed_voxels=[un_indexed_voxels;doubt_indexed_voxels];
end

% if only consider the already-indexed voxels with low completeness or belong to small grains
if correct_already_indexed_only==1
    un_indexed_voxels = find_doubt_indexed_voxels(DS_in,0.35,0);
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
                        Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,RotAxisOffset,pixelysize,pixelzsize,detysize,detzsize, ...
                        BeamStopY,BeamStopZ,RecVolumePixel,tomo_scale,VoxSize,simap_data_flag);
    grow_time=toc;
    fprintf('The find_neighbor_time and the grow_time take %0.2f s and %0.2f s, respectively\n',find_neighbor_time,grow_time)
    fprintf('Iter %d: reconstructed volume fraction: %.4f = %d / %d \n',i,length(find(DS_out.GrainId>0))/length(find(DS_out.Mask==1)), ...
        length(find(DS_out.GrainId>0)),length(find(DS_out.Mask==1)))
end

