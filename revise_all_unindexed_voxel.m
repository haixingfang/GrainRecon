% revise indexing for all the unindexed voxels (including indexed 1 or 2 voxels) by attempting indexed neighboring orientations
% Feb 14, 2022
function [DS_out]=revise_all_unindexed_voxel(DS_in,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,minComp,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,RecVolumePixel,simap_data_flag,revise_all,search_radius)
if nargin<35 || ~exist('search_radius','var')
    search_radius=20; % radius of search region to find candidate orientations [pixel]
    fprintf('search_radius = %d pixels\n',search_radius)
end
DS_out=DS_in;

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

% indices for un-indexed voxels
un_indexed_voxels=[];
[un_indexed_voxels(:,1), un_indexed_voxels(:,2), un_indexed_voxels(:,3)]=ind2sub(size(DS_in.VisitFlag), ...
        find(DS_in.GrainId==0 & DS_in.Mask==1));
if ~isempty(un_indexed_voxels)
    fprintf('find %d un-indexed voxels\n',length(un_indexed_voxels(:,1)))
else
    fprintf('find 0 un-indexed voxels\n')
end

% indices for indexed voxels yielding less than n-voxel size grains, n=2 by default
doubt_indexed_voxels = find_doubt_indexed_voxels(DS_in,minComp,1);

% all voxels waiting for revision
% revise_all=0;
if revise_all==1
    wait_revising_voxels=[doubt_indexed_voxels;un_indexed_voxels];
else
    wait_revising_voxels=doubt_indexed_voxels;
end
count=0;
if ~isempty(wait_revising_voxels)
    reindex_out=zeros(length(wait_revising_voxels(:,1)),5);
    sprintf('re-indexing and revising %d voxels (note: about 0.2 s for one voxel)...',length(wait_revising_voxels(:,1)))
    tic
    if length(wait_revising_voxels(:,1))>30000
        parfor i=1:length(wait_revising_voxels(:,1))
            id_vol=zeros(size(DS_in.GrainId));
            xn=wait_revising_voxels(i,1);
            yn=wait_revising_voxels(i,2);
            zn=wait_revising_voxels(i,3);
            id_vol(xn,yn,zn)=1;
            id_dismap=bwdist(id_vol);
            id_neigb=find(id_dismap>0 & id_dismap<search_radius); % distances follow: 1, sqrt(2), sqrt(3), 2, sqrt(5), ...
            id_neigb=unique(DS_in.GrainId(id_neigb));
            id_neigb=id_neigb(id_neigb>0); % neighbors ID
            pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*tomo_scale.VoxSize(1)+tomo_scale.Center'; % [mm]
            if simap_data_flag==1
                pos(1)=-pos(1);
                pos(2)=-pos(2);
            end
            [Comp, Dis_median,ind]=revise_all_unindexed_voxel_fun(id_neigb,DS_in,proj_bin_bw,pos,rot_angles, ...
                S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            if Comp(ind)>=DS_in.Completeness(xn,yn,zn)
                reindex_out(i,:)=[i,Comp(ind),Dis_median(ind),id_neigb(ind),ind];
            end
            fprintf('Voxel %d: old orientation (C = %.3f, grain_id = %d) vs new one (C = %.3f, grain_id = %d).\n', ...
                    i,DS_in.Completeness(xn,yn,zn),DS_in.GrainId(xn,yn,zn),Comp(ind),id_neigb(ind))
        end
    else
        for i=1:length(wait_revising_voxels(:,1))
            id_vol=zeros(size(DS_in.GrainId));
            xn=wait_revising_voxels(i,1);
            yn=wait_revising_voxels(i,2);
            zn=wait_revising_voxels(i,3);
            id_vol(xn,yn,zn)=1;
            id_dismap=bwdist(id_vol);
            id_neigb=find(id_dismap>0 & id_dismap<search_radius); % distances follow: 1, sqrt(2), sqrt(3), 2, sqrt(5), ...
            id_neigb=unique(DS_in.GrainId(id_neigb));
            id_neigb=id_neigb(id_neigb>0); % neighbors ID
            pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*tomo_scale.VoxSize(1)+tomo_scale.Center'; % [mm]
            if simap_data_flag==1
                pos(1)=-pos(1);
                pos(2)=-pos(2);
            end
            [Comp, Dis_median,ind]=revise_all_unindexed_voxel_fun(id_neigb,DS_in,proj_bin_bw,pos,rot_angles, ...
                S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            if Comp(ind)>=DS_in.Completeness(xn,yn,zn)
                reindex_out(i,:)=[i,Comp(ind),Dis_median(ind),id_neigb(ind),ind];
            end
            fprintf('Voxel %d: old orientation (C = %.3f, grain_id = %d) vs new one (C = %.3f, grain_id = %d).\n', ...
                    i,DS_in.Completeness(xn,yn,zn),DS_in.GrainId(xn,yn,zn),Comp(ind),id_neigb(ind))
        end
    end
    revise_time=toc;
    
    % update DS_out
    for i=1:length(wait_revising_voxels(:,1))
        xn=wait_revising_voxels(i,1);
        yn=wait_revising_voxels(i,2);
        zn=wait_revising_voxels(i,3);
        if reindex_out(i,2)>=DS_in.Completeness(xn,yn,zn)
            ind0=find(DS_in.GrainId==reindex_out(i,4));           
            DS_out.Completeness(xn,yn,zn)=reindex_out(i,2);
            DS_out.GrainId(xn,yn,zn)=reindex_out(i,4);
            DS_out.Rodrigues(:,xn,yn,zn)=DS_in.Rodrigues(:,ind0(1));
            DS_out.EulerZXZ(:,xn,yn,zn)=DS_in.EulerZXZ(:,ind0(1));
            DS_out.IPF001(:,xn,yn,zn)=DS_in.IPF001(:,ind0(1));
            DS_out.Dismedian(xn,yn,zn)=reindex_out(i,3);
            DS_out.VisitFlag(xn,yn,zn)=0.8; % not 1, meaning it is indexed by post processing
            count=count+1;
        end
    end
    fprintf('%d / %d voxels have been re-indexed and successfully revised: taking %0.2f s.\n', ...
        count,length(wait_revising_voxels(:,1)),revise_time)
else
    sprintf('No need to revise any voxels.')
end
end

function [Comp, Dis_median,ind]=revise_all_unindexed_voxel_fun(id_neigb,DS_in,proj_bin_bw,pos,rot_angles, ...
            S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

Comp=[];
Dis_median=[];
for j=1:size(id_neigb,1)
    UU=euler2u(DS_in.EulerAngle(id_neigb(j),1)*pi/180,DS_in.EulerAngle(id_neigb(j),2)*pi/180,DS_in.EulerAngle(id_neigb(j),3)*pi/180);
    [Nr_simu,Nr_intersect,dis_median]=forward_comp(UU,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
    Comp=[Comp;Nr_intersect/Nr_simu];
    Dis_median=[Dis_median;dis_median];
end
[~, ind] = max(Comp);
end
