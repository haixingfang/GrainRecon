% grain mapping engine with parallel computing orientations for for each seed
% step 1: calculate matched number of G vectors for each orientation for each seeding voxel position
% step 2: select orientation candidates based on the number of matched cases, the more the better
% step 3: calculation for each orientation contained in the local gridding space around each candidate 
% step 4: choose the maximum, and fit to maximize the completeness and find the optimal orientation
% step 5: fit the grain center of mass according to the spots paired by the forward simulation
% step 6: grow regions around the seeding voxel by controlling the relative decrease of the completeness
% compared to GrainMapping_engine_v2.m, it adds more choices of values for VisitFlag
% when growing voxels, also compare the number of intersections because wrong indexing yielding the same completeness value is possible
% April 8, 2022
function [indexing_final,pos_list_new,DS_out,Nr_indexed_voxel]=GrainMapping_engine_gpu_comp(pos_list,OR,RotDet, ...
        Spots,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
        Lsam2sou,Lsam2det,TrustComp,minComp,minEucDis,dety00,detz00,P0y,P0z, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        drop_off,maxD,DS,tomo_scale,VoxSize,RecVolumePixel,FirstGrainID,simap_data_flag, ...
        maxDmedian,Nr_seed,OutputFolder,iter,check_fit_all)
% % for testing
% pos_list=pos_seed;
% if length(FirstGrainID)>1
%     FirstGrainID=FirstGrainID(iter);
% end
% % Spots=SpotsForIndex;
% check_fit_all=0;

% initialzie
PhaseId=DS.PhaseId;
Mask=DS.Mask;
Completeness=DS.Completeness;
GrainId=DS.GrainId;
Rodrigues=DS.Rodrigues;
IPF001=DS.IPF001;
EulerZXZ=DS.EulerZXZ;
Dismedian=DS.Dismedian;
Ninter=DS.Ninter;
Icorr=DS.Icorr;
VisitFlag=DS.VisitFlag;
Nr_indexed_voxel=zeros(1,length(pos_list(:,1)));
indexing_final0=zeros(length(pos_list(:,1)),18);
indexing_final=zeros(length(pos_list(:,1)),18);
confident_index=zeros(1,length(pos_list(:,1)));
replace_center_count=zeros(length(pos_list(:,1)),1);
pos_list_new=zeros(length(pos_list(:,1)),3);
SeedIndexing=zeros(length(pos_list(:,1)),15);
candidate_id=zeros(length(pos_list(:,1)),1);
fit_all_flag_track=zeros(length(pos_list(:,1)),1);
if check_fit_all==1
    ind_GrainId=[];
    [ind_GrainId(:,1),ind_GrainId(:,2),ind_GrainId(:,3)]=ind2sub(size(DS.GrainId),find(DS.GrainId>0));
end
for i=1:length(pos_list(:,1))
    stop_regional_recon=0;
    pos_indexing=pos_list(i,:); % position for indexing [mm]
    sprintf('Iter %d: indexing processing %d / %d, voxel position [%.4f %.4f %.4f]',iter,i,length(pos_list(:,1)),pos_indexing)
    while stop_regional_recon~=1 && replace_center_count(i)<=3
        pos_indices=pos2ind(i,pos_indexing,tomo_scale,VoxSize,RecVolumePixel,size(Completeness),simap_data_flag);
        if replace_center_count(i)~=0
            fprintf('Iter %d: indexing processing %d / %d, updated voxel position [%.4f %.4f %.4f].\n',iter,i,length(pos_list(:,1)),pos_indexing)
        end
        if (Completeness(pos_indices(1),pos_indices(2),pos_indices(3))<TrustComp || (replace_center_count(i)>0 && ...
                Completeness(pos_indices(1),pos_indices(2),pos_indices(3))>=TrustComp)) && ...
                VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))~=1
%         sprintf('Indexing processing %d / %d, voxel position [%.3f %.3f %.3f]',i,length(pos_list(:,1)),pos_indexing)
        tic
        if replace_center_count(i)==0
%             Output = ForwardCalc_OR_space_Gt_cuda(pos_indexing,OR,S,B,Ahkl,nrhkl,hkl_family_square,d_possible, ...
%                         Spots,rot_angles,RotDet,rot_start,rot_step,Lsam2sou,Lsam2det, ...
%                         P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize);
            Output = ForwardCalc_OR_space_comp_cuda(pos_indexing,OR,S,B,Ahkl,nrhkl, ...
                        proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                        P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                        BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max);
            if sum(all(Output(:,2:4)==0))==3
                error('error: cuda_forward_comp.cu does not return correct results !');
            end
%             [Output,~]=ForwardCalc_OR_space_Gt(OR.q_b,OR.q_c,OR.q_d,size(OR.q_b,1),RotDet, ...
%                 Spots,pos_indexing,rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,hkl_family, ...
%                 hkl_family_square,d_possible,Glen_possible,Lsam2sou,Lsam2det,dety00,detz00, ...
%                 P0y,P0z,pixelysize,pixelzsize,dety0,detz0,thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ);
        else
            % 8: 730 local orientations
            % 4: 126 local orientations
            % 3: 65 local orientations
            OR_local=ori_local_sampling(indexing_final_temp(15:18),4);
            Output = ForwardCalc_OR_space_comp_cuda(pos_indexing,OR_local,S,B,Ahkl,nrhkl, ...
                        proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                        P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                        BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max);
%             Output = ForwardCalc_OR_space_Gt_cuda(pos_indexing,OR_local,S,B,Ahkl,nrhkl,hkl_family_square,d_possible, ...
%                         Spots,rot_angles,RotDet,rot_start,rot_step,Lsam2sou,Lsam2det, ...
%                         P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize);
%             [Output,~]=ForwardCalc_OR_space_Gt(OR_local.q_b,OR_local.q_c,OR_local.q_d,size(OR_local.q_b,1),RotDet, ...
%                 Spots,pos_indexing,rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,hkl_family, ...
%                 hkl_family_square,d_possible,Glen_possible,Lsam2sou,Lsam2det,dety00,detz00, ...
%                 P0y,P0z,pixelysize,pixelzsize,dety0,detz0,thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ);
        end
        OR_time=toc;
        index_candidates=double(Output);
%         [index_candidates,~]=sortrows(index_candidates,2,'descend'); % number of matched G vectors
        [index_candidates,~]=sortrows(index_candidates,3,'descend'); % completeness, 
%         [index_candidates,~]=sortrows(index_candidates,4,'ascend'); % median distance

        if length(index_candidates(:,1))<Nr_seed
            select_index3=1:length(index_candidates(:,1));
        else
            select_index3=1:Nr_seed;
        end
        index_seeds=index_candidates(select_index3,:);
        
        % fit to obtain orientation by maximizing the completeness value,
        % constrain fitting
        tic
        Spots_pair{i}=zeros(length(index_seeds(:,1)),11);
%         indexing_fit{i}=zeros(length(index_seeds(:,1)),18);
        if replace_center_count(i)==0
            if check_fit_all==1 && min(sqrt((ind_GrainId(:,1)-pos_indices(1)).^2+(ind_GrainId(:,2)-pos_indices(2)).^2+ ...
                    (ind_GrainId(:,3)-pos_indices(3)).^2))>5
                fit_all_flag=1;
                fprintf('Fit all OR candidates for seed %d.\n',i)
            else
                fit_all_flag=0;
            end
%             [indexing_final(i,:),indexing_fit{i},confident_index(i),Spots_pair{i}]=fit_OR_constrain_with_local_comp(i,index_seeds,RotDet, ...
%                 proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
%                 thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
%                 pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            [indexing_final(i,:),confident_index(i),Spots_pair{i}]=fit_OR_constrain_with_local_comp_cuda(i,index_seeds,RotDet, ...
                proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
                thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ,fit_all_flag);
            indexing_final0(i,:)=indexing_final(i,:);
            fprintf('Candidate %d is identified for yielding the best solution.\n',confident_index(i))
			candidate_id(i)=confident_index(i);
            fit_all_flag_track(i)=fit_all_flag;
            
            if ~isempty(Spots_pair{i})
			    % fit the grain center of mass based on centers of the paired experimental spots
			    [FitOutput,fval]=fit_grain_COM(Spots_pair{i},S,B,Lsam2sou,Lsam2det,P0y,P0z, ...
				    pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet);
			    if sqrt(sum((FitOutput(1:3)-indexing_final0(i,2:4)).^2))/VoxSize>maxD % update the voxel position, being closer to the grain COM
				    fprintf('The spot centroid difference is %.2f pixels after fitting the grain COM.\n',fval)
                    pos_indices_temp=pos2ind(i,FitOutput(1:3),tomo_scale,VoxSize,RecVolumePixel,size(Completeness),simap_data_flag);
                    if Completeness(pos_indices_temp(1),pos_indices_temp(2),pos_indices_temp(3))<indexing_final(i,7) ...
                            && Mask(pos_indices_temp(1),pos_indices_temp(2),pos_indices_temp(3))==1
                        indexing_final(i,2:4)=FitOutput(1:3);
                        if length(FitOutput)>3
                            indexing_final(i,9:11)=FitOutput(4:6);
                            indexing_final(i,12:14)=angle2rod(FitOutput(4)*pi/180,FitOutput(5)*pi/180,FitOutput(6)*pi/180,'ZXZ'); % Rodrigues
                            indexing_final(i,15:18)=rod2quat(indexing_final(i,12:14)); % quaternion
                        end
                        UU=euler2u(indexing_final(i,9)*pi/180,indexing_final(i,10)*pi/180,indexing_final(i,11)*pi/180);
                        [indexing_final(i,5),indexing_final(i,6),indexing_final(i,8)]=forward_comp(UU,proj_bin_bw,indexing_final(i,2:4),rot_angles,S,B,Ahkl,nrhkl, ...
                            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
                        indexing_final(i,7)=indexing_final(i,6)/indexing_final(i,5);
                        fprintf('Update the voxel position from [%.4f %.4f %.4f] to [%.4f %.4f %.4f] mm,\n euler angles from [%.2f %.2f %.2f] to [%.2f %.2f %.2f] degrees,\n completeness from %.3f to %.3f.\n', ...
                            pos_indexing,indexing_final(i,2:4),indexing_final0(i,9:11),indexing_final(i,9:11),indexing_final0(i,7),indexing_final(i,7))
                        pos_indexing=FitOutput(1:3);
                        pos_indices=pos2ind(i,pos_indexing,tomo_scale,VoxSize,RecVolumePixel,size(Completeness),simap_data_flag);
                    end
                end
            end
        else
%             [indexing_final(i,:),indexing_fit{i},confident_index(i),Spots_pair{i}]=fit_OR_constrain_local(i,index_seeds,RotDet, ...
%                 proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
%                 Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
%                 pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            [indexing_final(i,:),confident_index(i),Spots_pair{i}]=fit_OR_constrain_with_local_comp_cuda(i,index_seeds,RotDet, ...
                proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
                thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        end
        fit_time=toc;
        fprintf('The OR calc and the fitting take %0.2f s and %0.2f s, respectively.\n',OR_time,fit_time)
        fprintf('Indexed completeness is %.3f and Euler angles are [%.2f, %.2f, %.2f] degrees.\n',indexing_final(i,[7 9:11]))
        
        if indexing_final(i,7)<minComp || indexing_final(i,7)<Completeness(pos_indices(1),pos_indices(2),pos_indices(3)) ...
                || indexing_final(i,8)>maxDmedian
            stop_regional_recon=1;
            pos_list_new(i,:)=pos_indexing;
            VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))=0.8; % it means this voxel is a failed seeding voxel
            fprintf('Not to grow because the completeness (%.3f) is too low (< %.3f or < %.3f) or the median D (%.2f) is too large (>%.1f pixel).\n', ...
                indexing_final(i,7),minComp,Completeness(pos_indices(1),pos_indices(2),pos_indices(3)),indexing_final(i,8),maxDmedian)
        else
            % update DS and grow around this voxel if its completeness value is
            % higher than the old one
            fprintf('Grow the region because the completeness (%.3f) is >= %.3f and the median D (%.2f) is <= %.1f pixels.\n', ...
                indexing_final(i,7),minComp,indexing_final(i,8),maxDmedian)
            indexing_final_temp=indexing_final(i,:);
            center=pos_indices;
            replace_center=0;
            if indexing_final_temp(7)>=Completeness(pos_indices(1),pos_indices(2),pos_indices(3))
                Completeness(pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final_temp(7);
                indexed_comp=indexing_final_temp(7);
                % grow indexed region
                tic
                UU=quaternion2U(indexing_final_temp(15:18));
                if ~isempty(Spots_pair{i})
                    Spot_L = mean(Spots_pair{i}(:,21))+2*std(Spots_pair{i}(:,21)); % bounding area of the spot [pixel^2]
                    Spot_L = sqrt(Spot_L/pi); % bounding half length of the spot [pixel]
                else
                    Spot_L = 20;
                end
                if Spot_L <= 2
                    fprintf("Spot_L = %.2f pixels, grow by connecting voxels (cpu computation).\n", Spot_L)
                    [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region(Completeness,Mask,Dismedian,Ninter,pos_indices,drop_off,maxD, ...
                        indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                            RotDet,thetamax,lambda_min,lambda_max, ...
                        Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian);
                else
                    fprintf("Spot_L = %.2f pixels, grow within a bounding volume (gpu-cuda computation).\n", Spot_L)
%                     [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region_parallel(Completeness,Mask,Dismedian,Ninter,pos_indices,Spot_L,drop_off,maxD, ...
%                         indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
%                             RotDet,thetamax,lambda_min,lambda_max, ...
%                         Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
%                         pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
%                         RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian);
                    [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region_cuda(Completeness,Mask,Dismedian,Ninter,pos_indices,Spot_L,drop_off,maxD, ...
                        indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                            RotDet,thetamax,lambda_min,lambda_max, ...
                        Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                        pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian);
                end
                ind=sub2ind(size(J),indexed_indices(:,1),indexed_indices(:,2),indexed_indices(:,3));
                Nr_indexed_voxel(i)=length(ind);
    %             PhaseId(ind)=1;
                GrainId(ind)=i+FirstGrainID-1;
    %             [ind1 ind2 ind3]=ind2sub(size(J),find(PhaseId>0));
                Rodrigues(:,ind)=repmat(indexing_final_temp(12:14)',1,length(ind));
                IPF001(:,ind)=repmat(indexing_final_temp(9:11)'./norm(indexing_final_temp(9:11)),1,length(ind)); % colored by Euler angles
                EulerZXZ(:,ind)=repmat(indexing_final_temp(9:11)',1,length(ind));
                VisitFlag(ind)=0.5;
                VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))=1;
                Completeness(ind)=Completeness_out(ind);
                Dismedian(ind)=DisMedian_out(ind);
                Ninter(ind)=Ninter_out(ind);
                grow_time=toc;
                fprintf('%d voxels have grown around the seed %d and the growth takes %0.2f s.\n',length(ind)-1,i,grow_time)
                fprintf('center updated from [%d %d %d] to [%d %d %d].\n',pos_indices,center)
                if replace_center==1 && GrainId(center(1),center(2),center(3))~=i+FirstGrainID-1 % exclude the possibility to find a new center belong to other grains
%                     replace_center=0;
%                     center=pos_indices;
                    AllID=[];
                    [AllID(:,1),AllID(:,2),AllID(:,3)]=ind2sub(size(GrainId),find(GrainId==i+FirstGrainID-1));
                    [~,ind_ALLID]=min(sqrt(sum((AllID-center).^2,2)));
                    center=AllID(ind_ALLID,:); % find the voxel which is closest to the true center
%                     replace_center=0;
%                     sprintf('center will not be re-updated from [%d %d %d] to [%d %d %d]',pos_indices,center)
                    if sqrt(sum((center-pos_indices).^2))>maxD
                        fprintf('center re-updated from [%d %d %d] to [%d %d %d].\n',pos_indices,center)
                        replace_center=1;
                    else
                        fprintf('center [%d %d %d] is close enough to previous [%d %d %d], stop center updating !\n',pos_indices,center)
                        replace_center=0;
                    end
                end
            end

            if replace_center==0
                stop_regional_recon=1;
                pos_list_new(i,:)=pos_indexing;           
            else
                % go back to initial state
                replace_center_count(i)=replace_center_count(i)+1;
%                 Nr_indexed_voxel(i)=0;
%                 GrainId(ind)=0;
%                 Rodrigues(:,ind)=repmat([0 0 0]',1,length(ind));
%                 IPF001(:,ind)=repmat([0 0 0]',1,length(ind));
%                 EulerZXZ(:,ind)=repmat([0 0 0]',1,length(ind));
%                 VisitFlag(ind)=0.1; % This voxel had grown before but has to go back to 'unindexed' because of the change of the weighted center
%                 Completeness(ind)=0;
%                 Dismedian(ind)=10;
%                 Ninter(ind)=0;
                pos_list_new(i,:)=((center+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
                if simap_data_flag==1
                    pos_list_new(i,1)=-pos_list_new(i,1);
                    pos_list_new(i,2)=-pos_list_new(i,2);
                end
                pos_indexing=pos_list_new(i,:);
            end
        end
        else
            stop_regional_recon=1;
            pos_list_new(i,:)=pos_indexing; 
        end
    end
    if i>1
        SeedIndexing(i,:)=[indexing_final(i,1:11) Nr_indexed_voxel(i) ...
            length(find(GrainId>0))-length(find(DS_out.GrainId>0)) candidate_id(i) fit_all_flag_track(i)];
    else
        SeedIndexing(i,:)=[indexing_final(i,1:11) Nr_indexed_voxel(i) ...
            length(find(GrainId>0))-length(find(DS.GrainId>0)) candidate_id(i) fit_all_flag_track(i)];
    end
    DS_out.PhaseId=PhaseId;
    DS_out.Mask=Mask;
    DS_out.Completeness=Completeness;
    DS_out.GrainId=GrainId;
    DS_out.Rodrigues=Rodrigues;
    DS_out.EulerZXZ=EulerZXZ;
    DS_out.IPF001=IPF001;
    DS_out.Dismedian=Dismedian;
    DS_out.Ninter=Ninter;
    DS_out.Icorr=Icorr;
    DS_out.VisitFlag=VisitFlag;
    if mod(i,50)==0 || i==length(pos_list(:,1))
        save(fullfile(OutputFolder,strcat(num2str(iter), 'DS.mat')),'DS_out','-v7.3');
        save(fullfile(OutputFolder,[strcat('SeedingIndexing',num2str(iter)) '.mat']),'SeedIndexing');
    end
    sprintf('Reconstructed volume fraction: %.4f = %d / %d ',length(find(DS_out.GrainId>0))/length(find(DS_out.Mask==1)), ...
        length(find(DS_out.GrainId>0)),length(find(DS_out.Mask==1)))
    % add on June 10, 2022
    if check_fit_all~=1 && i>200 && ~isempty(find(SeedIndexing(i-180:i,1)>0))
        if length(find(SeedIndexing(i-180:i,13)>0))/length(find(SeedIndexing(i-180:i,1)>0))<0.05
            check_fit_all=1;
            ind_GrainId=[];
            [ind_GrainId(:,1),ind_GrainId(:,2),ind_GrainId(:,3)]=ind2sub(size(DS.GrainId),find(DS.GrainId>0));
            sprintf("Grain mapping engine will check the need of fit_all mode for seeds.")
        end
    end
    if check_fit_all==1 && i>200 && ~isempty(find(SeedIndexing(i-50:i,1)>0))
        if length(find(SeedIndexing(i-50:i,13)>0))/length(find(SeedIndexing(i-50:i,1)>0))<0.05
            break;
        end
    end
end
