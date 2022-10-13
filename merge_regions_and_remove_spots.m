function [Spots_keep,DS_merge]=merge_regions_and_remove_spots(DS,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber, ...
                hkl_square,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                tomo_scale,VoxSize,RecVolumePixel,simap_data_flag,dis_tolerant,min_misori)
% merge regions that have smaller misorientation than pre-defined threshold value
% Feb 9, 2022
mtex_avail=0;  % availability of mtex toolbox, 1: yes; 0: no (default).
if mtex_avail~=0
    cs = crystalSymmetry('cubic'); % 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'trigonal'
else
    cs = [];
end
if nargin<37
    dis_tolerant = 10; % tolorance spot distance to accept 'assigned spots' [pixel]
    min_misori = 0.5; % recommend to be 0.5 [deg]
end
if nargin<38
    min_misori = 0.5; % recommend to be 0.5 [deg]
end
[DS_merge,Ngrain,~,~,~]=merge_and_update_completeness(DS,mtex_avail, ...
    cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
iter_merge=1;
Ngrain_iter(iter_merge)=Ngrain;
iter_flag=1;  % iterative merging
if iter_flag==1
    stop_iter=0;
    while stop_iter~=1
        [DS_merge,Ngrain,~,~,~]=merge_and_update_completeness(DS_merge,mtex_avail, ...
            cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
        iter_merge=iter_merge+1;
        Ngrain_iter(iter_merge)=Ngrain;
        if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=0
            stop_iter=1;
        end
    end
end

% revise indexing for the unindexed single voxel
[DS_merge]=revise_single_unindexed_voxel(DS_merge);
DS_merge.Completeness(DS_merge.Completeness>1)=1;

DS_merge.SeedID = 1:max(DS_merge.GrainId(:));
for i = DS_merge.SeedID
    [x,y,z] = ind2sub(size(DS_merge.GrainId),find(DS_merge.GrainId == i));
    X = mean(x);
    Y = mean(y);
    Z = mean(z);
    DS_merge.Coord(i,:) = [X,Y,Z]; % each row represents coodinate for each grain
    DS_merge.nVox(i,1) = length(find(DS_merge.GrainId == i)); % number of voxels for each grain
    DS_merge.EulerAngle(i,:)=DS_merge.EulerZXZ(:,x(1),y(1),z(1))';
end

% find spots corresponding to indexed grains
remove_spotID=[];
for i=1:length(DS_merge.SeedID)
    if DS_merge.nVox(i)>0
        pos=((DS_merge.Coord(i,:)+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
        if simap_data_flag==1
            pos(1:2)=-pos(1:2);
        end
        U=euler2u(DS_merge.EulerAngle(i,1)*pi/180,DS_merge.EulerAngle(i,2)*pi/180,DS_merge.EulerAngle(i,3)*pi/180);
        [~,~,~,~,HittedSpots]=index_verify_v3(U,proj_bin_bw,Spots,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
	    pixelnr_pred=(3*DS_merge.nVox(i)/(4*pi))^(1/3)*2*VoxSize./(mean([pixelysize pixelzsize]))*Lsam2det/Lsam2sou; % [pixel]
        if ~isempty(HittedSpots)
            for j=1:length(HittedSpots(:,1))
                ProjectNo=(HittedSpots(j,17)-rot_start)/rot_step+1; % project no.
                SpotID=find(Spots{ProjectNo}(:,9)>(HittedSpots(j,18)-0.1) & Spots{ProjectNo}(:,9)<(HittedSpots(j,18)+0.1) ...
                    & Spots{ProjectNo}(:,10)>(HittedSpots(j,19)-0.1) & Spots{ProjectNo}(:,10)<(HittedSpots(j,19)+0.1));
                if ~isempty(SpotID) && HittedSpots(j,20)<dis_tolerant
		    if abs(pixelnr_pred-sqrt(4*Spots{ProjectNo}(SpotID,11)/pi))/sqrt(4*Spots{ProjectNo}(SpotID,11)/pi)<1/3
                    	remove_spotID=[remove_spotID;ProjectNo SpotID HittedSpots(j,20)];
		    end
                end
            end
        end
    end
    i;
end
[~,ia,~] = unique(remove_spotID(:,1:2),'rows');
remove_spotID=remove_spotID(ia,:);

%%%%%%%%%%%%%%%%%% remove the assigned spots from the spot list
for i=1:length(Spots)
    to_remove_spotID=remove_spotID(remove_spotID(:,1)==i,:);
    if ~isempty(to_remove_spotID)
        keep_spotID=setdiff([1:length(Spots{i}(:,1))]',to_remove_spotID(:,2),'rows');
        Spots_keep{i}=Spots{i}(keep_spotID,:);
    else
        Spots_keep{i}=Spots{i};
    end
    if i==1 || i==length(Spots) 
    	sprintf('%d / %d spots are removed from the spot list for projection %d', ...
             length(to_remove_spotID(:,1)),length(Spots{i}(:,1)),i)
    end
end
