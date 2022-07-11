% define grains: merge regions that have smaller misorientation than pre-defined threshold value
% input: DS, output from the previous index_grow process
%        mtex_avail = 1 or 0
%        cs: crystalSymmetry: 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'trigonal'
%        min_misori = minimum misorientation to identify grains, 0.5 deg by default 
% output: DS_merge; i - number of identified grains; id0 - number of regions
% Using mtex toolbox for calculating the misorientation angle, because function <misori2> might not be capable to cover certain symmetry
function [DS_merge,i,id0,Inherit_region_nr,CentroidComp]=merge_and_update_completeness_mtex(DS,mtex_avail,cs,min_misori,proj_bin_bw, ...
        Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        tomo_scale,VoxSize,RecVolumePixel,simap_data_flag)
    
% % for testing
% mtex_avail=1;  % availability of mtex toolbox, 1: yes; 0: no (default).
% if mtex_avail~=0
%     cs = crystalSymmetry('cubic'); % 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'trigonal'
% end
% min_misori = 0.5; % [deg]

if isempty(mtex_avail)
    mtex_avail=0;
    cs=[];
else
    cs = crystalSymmetry('cubic');
end
if isempty(min_misori)
    min_misori=0.5;
end
DS_merge=DS;
DS_merge.GrainId=zeros(size(DS.GrainId));
DS_merge.Rodrigues=zeros(size(DS.Rodrigues));
DS_merge.EulerZXZ=zeros(size(DS.EulerZXZ));
DS_merge.IPF001=zeros(size(DS.IPF001));

id=unique(reshape(DS.GrainId,[],1));
id=id(id>0);
id0=id;
stop_merge=0;
i=0;
update_completeness=1; % 1: to update completeness; 0: not to update completeness
while stop_merge~=1
    i=i+1;
    id_vol=DS.GrainId==id(1); % must be one single connected region
    id_dismap=bwdist(id_vol);
    id_neigb=find(id_dismap>0 & id_dismap<2); % distances follow: 1, sqrt(2), sqrt(3), 2, sqrt(5), ...
    id_neigb=unique(DS.GrainId(id_neigb));
    id_neigb=id_neigb(id_neigb>0); % neighbors ID
    id_indices=[];
    [id_indices(:,1), id_indices(:,2), id_indices(:,3)]=ind2sub(size(id_vol),find(id_vol==1));
    r0=DS.Rodrigues(:,id_indices(1,1),id_indices(1,2),id_indices(1,3))';
%     if all(r0==0)
%         r0=DS.Rodrigues(:,id_indices(round(length(id_indices(:,1))/2),1), ...
%             id_indices(round(length(id_indices(:,1))/2),2),id_indices(round(length(id_indices(:,1))/2),3))';
%     end
    q0=rod2quat(r0);
    U0=quaternion2U(rod2quat(r0));
    euler_angles0=u2euler_corr(U0);
    V0=length(find(DS.GrainId==id(1))); % number of voxels
    [~,Cmax0]=max(DS.Completeness(sub2ind(size(DS.GrainId),id_indices(:,1),id_indices(:,2),id_indices(:,3))));
    Cmax0=id_indices(Cmax0,:);          % indices for the maximum completeness
    mergeInd=id_indices;                % indices for merging regions
    if ~isempty(id_neigb)
        clear r1 q1 V1 ang ax euler_angles1 Cmax1 update_completeness_factor;
        for j=1:length(id_neigb)
           neigb_indices{j}=[];
           [neigb_indices{j}(:,1), neigb_indices{j}(:,2), neigb_indices{j}(:,3)]=ind2sub(size(DS.GrainId),find(DS.GrainId==id_neigb(j)));
           r1(j,:)=DS.Rodrigues(:,neigb_indices{j}(1,1),neigb_indices{j}(1,2),neigb_indices{j}(1,3))';
           q1(j,:)=rod2quat(r1(j,:));
           V1(j)=length(find(DS.GrainId==id_neigb(j))); % number of voxels
           U1=quaternion2U(rod2quat(r1(j,:)));
           euler_angles1(j,:)=u2euler_corr(U1);
%            [ang(j), ax(j,:)] = misori2(U0,U1);
           
           % using mtex, Feb 11, 2022
           rot0 = rotation('Euler',euler_angles0(1)*degree,euler_angles0(2)*degree,euler_angles0(3)*degree);
           rot1 = rotation('Euler',euler_angles1(j,1)*degree,euler_angles1(j,2)*degree,euler_angles1(j,3)*degree);
           o0=orientation(rot0,cs);
           o1=orientation(rot1,cs);
           ang(j)=angle(o0,o1)/degree;
        end
        ID_for_merge=id_neigb(find(ang<=min_misori));   % neighboring ID for merging
        ID_for_merge_indices=find(ang<=min_misori);
        if ~isempty(ID_for_merge)
            Inherit_region_nr{i} = [id(1) ID_for_merge'];
            weights=[V0 V1(ID_for_merge_indices)]./sum([V0 V1(ID_for_merge_indices)]); % weights for updating the orientation and completeness
            if mtex_avail==1
                % using mtex
                rot = rotation('Euler',[euler_angles0(1) euler_angles1(ID_for_merge_indices,1)']'*degree, ...
                    [euler_angles0(2) euler_angles1(ID_for_merge_indices,2)']'*degree, ...
                    [euler_angles0(3) euler_angles1(ID_for_merge_indices,3)']'*degree);
                ori=orientation(rot,cs);
        %         [m, q, lambda, V] = mean(ori,'weights',weights)
                [EulerZXZ_new, ~, ~, ~] = mean(ori,'weights',weights);
                EulerZXZ_new = [EulerZXZ_new.phi1 EulerZXZ_new.Phi EulerZXZ_new.phi2]*180/pi; % [deg]
                [q_new,r_new]=ang2quat(EulerZXZ_new);
            else
                % self defined:
                % 1) using weighted average of Rodrigues => 0.013 deg difference with mtex
                % 2) using weighted average of Quaternions => 0.06 deg difference with mtex
                r_new=([r0;r1(ID_for_merge_indices,:)]'*weights')';       % update Rodrigues vector
                q_new=rod2quat(r_new); % quaternion
    %             q_new=weights*[q0;q1(ID_for_merge_indices,:)];
                q_new=q_new/norm(q_new);
                EulerZXZ_new=u2euler_corr(quaternion2U(q_new));
    %             misori2(euler2u(EulerZXZ_new(1)*pi/180,EulerZXZ_new(2)*pi/180,EulerZXZ_new(3)*pi/180), ...
    %                 quaternion2U(q_new))
            end
            for k=1:length(ID_for_merge)
                mergeInd=[mergeInd;neigb_indices{ID_for_merge_indices(k)}]; % indices for merging regions
                [~,Cmax1_ind]=max(DS.Completeness(sub2ind(size(DS.GrainId), ...
                    neigb_indices{ID_for_merge_indices(k)}(:,1),neigb_indices{ID_for_merge_indices(k)}(:,2), ...
                    neigb_indices{ID_for_merge_indices(k)}(:,3))));
                Cmax1(k,:)=neigb_indices{ID_for_merge_indices(k)}(Cmax1_ind,:);
            end
        else
            Inherit_region_nr{i} = id(1);
            EulerZXZ_new=euler_angles0;
            r_new=r0;
        end
    
        % update DS_merge
        ind=sub2ind(size(DS.GrainId),mergeInd(:,1),mergeInd(:,2),mergeInd(:,3)); % linear indices
        DS_merge.GrainId(ind)=i;
        DS_merge.Rodrigues(:,ind)=repmat(r_new',1,length(ind));
        DS_merge.IPF001(:,ind)=repmat(EulerZXZ_new'./norm(EulerZXZ_new),1,length(ind)); % colored by Euler angles
        DS_merge.EulerZXZ(:,ind)=repmat(EulerZXZ_new',1,length(ind));
        U_new=quaternion2U(rod2quat(r_new));
        if update_completeness && ~isempty(ID_for_merge)
            % consider the weights
            for k=1:length(Cmax1(:,1))+1
                if k==1
                    pos_indices=Cmax0;
                else
                    pos_indices=Cmax1(k-1,:);
                end
                if simap_data_flag==1
                    pos_indexing=((pos_indices+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center';
                    pos_indexing(1)=-pos_indexing(1);
                    pos_indexing(2)=-pos_indexing(2);
                else
                    pos_indexing=((pos_indices+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center';
                end
                % indexing verification
                [Nr_simu,Nr_intersect,~]=forward_comp(U_new,proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
                    RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);      
                update_completeness_factor(k)=(Nr_intersect/Nr_simu)/DS.Completeness(pos_indices(1),pos_indices(2),pos_indices(3));
            end
            DS_merge.Completeness(ind)=DS.Completeness(ind)*sum(update_completeness_factor.*weights); % times the weighted factor
            if length(find(DS_merge.Completeness==Inf))>0
                error('Unexpected Inf value for the completeness is found.');
            elseif DS_merge.Completeness(ind)>1
                DS_merge.Completeness(ind)=1;
            end
        end
                   
        % remove the id which have been processed
        id=setdiff(id,[id(1);ID_for_merge]);
%         if ~isempty(ID_for_merge)
%             sprintf('%d grains identified. %d regions have been merged and %d regions waiting for merging ...',i,length(ID_for_merge)+1,length(id))
%         else
%             sprintf('%d grains identified. 0 region has been merged and %d regions waiting for merging ...',i,length(id))
%         end
    else
        % update DS_merge
        Inherit_region_nr{i} = id(1);
        U_new=quaternion2U(rod2quat(r0));
        ind=sub2ind(size(DS.GrainId),mergeInd(:,1),mergeInd(:,2),mergeInd(:,3)); % linear indices
        DS_merge.GrainId(ind)=i;
        DS_merge.Rodrigues(:,ind)=repmat(r0',1,length(ind));
        DS_merge.IPF001(:,ind)=repmat(euler_angles0'./norm(euler_angles0),1,length(ind)); % colored by Euler angles
        DS_merge.EulerZXZ(:,ind)=repmat(euler_angles0',1,length(ind));
        DS_merge.Completeness(ind)=DS.Completeness(ind);
        id=setdiff(id,id(1));
%         sprintf('%d grains identified. 0 region has been merged and %d regions waiting for merging ...',i,length(id))
    end
    
    % find the hitted spots
    pos_indices=round(median(mergeInd,1));
    AllID=[];
    [AllID(:,1),AllID(:,2),AllID(:,3)]=ind2sub(size(DS_merge.GrainId),ind);
    [~,ind0]=min(sqrt(sum((AllID-pos_indices).^2,2)));
    pos_indices=AllID(ind0,:); % find the voxel which is closest to the true center
    
    if simap_data_flag==1
        pos_indexing=((pos_indices+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center';
        pos_indexing(1)=-pos_indexing(1);
        pos_indexing(2)=-pos_indexing(2);
    else
        pos_indexing=((pos_indices+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center';
    end
%     % indexing verification
    [Nr_simu,Nr_intersect,~]=forward_comp(U_new,proj_bin_bw,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
    CentroidComp(i)=Nr_intersect/Nr_simu;
    
    % update the completeness again, Dec 9, 2021
    if DS_merge.Completeness(pos_indices(1),pos_indices(2),pos_indices(3))>0
        DS_merge.Completeness(ind)=DS_merge.Completeness(ind)*(CentroidComp(i)/DS_merge.Completeness(pos_indices(1),pos_indices(2),pos_indices(3)));
    end
    
    if isempty(id)
        stop_merge=1;
    end
end
sprintf('%d grains identified out of %d regions.',i,length(id0))

