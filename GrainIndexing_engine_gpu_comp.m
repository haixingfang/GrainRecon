% gpu_cuda computing only
% grain indexing engine with parallel computing orientations for for each seed
% step 1: calculate completeness and dis_median for each orientation for each seeding voxel position
% step 2: select orientation candidates based on the completeness, the higher the better
% step 3: calculation for each orientation contained in the local gridding space around each candidate 
% step 4: choose the maximum, and fit to maximize the completeness and find the optimal orientation
% April 19, 2022
function [indexing_final,DS_out,pos_list_new]=GrainIndexing_engine_gpu_comp(pos_list,OR,RotDet, ...
                    Spots,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
                    Lsam2sou,Lsam2det,TrustComp,minComp,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                    DS,tomo_scale,VoxSize,RecVolumePixel,FirstGrainID,simap_data_flag,maxDmedian, ...
                    maxD,Nr_seed,OutputFolder,iter)
% % for testing
% pos_list=pos_seed;
% if length(FirstGrainID)>1
%     FirstGrainID=FirstGrainID(iter);
% end
% Spots=SpotsForIndex;

% initialzie
pos_list_new=pos_list;
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

indexing_final0=zeros(length(pos_list(:,1)),18);
indexing_final=zeros(length(pos_list(:,1)),18);
confident_index=zeros(1,length(pos_list(:,1)));

fprintf('Progress in indexing ... ...\n');
for i=1:length(pos_list(:,1))
    pos_indexing=pos_list(i,:); % position for indexing [mm]
    pos_indices=pos2ind(i,pos_indexing,tomo_scale,VoxSize,RecVolumePixel,size(Completeness),simap_data_flag);
    sprintf('Iter %d: indexing processing %d / %d, voxel position [%.4f %.4f %.4f]',iter,i,length(pos_list(:,1)),pos_indexing)
    if Completeness(pos_indices(1),pos_indices(2),pos_indices(3))<TrustComp && VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))~=1
        tic
        Output = ForwardCalc_OR_space_comp_cuda(pos_indexing,OR,S,B,Ahkl,nrhkl, ...
                    proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                    P0y,P0z,RotAxisOffset,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                    BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max);
        if sum(all(Output(:,2:4)==0))==3
            error('error: cuda_forward_Gt.cu does not return correct results !');
        end
        OR_time=toc;
        index_candidates=double(Output);
        [index_candidates,~]=sortrows(index_candidates,3,'descend'); % completeness, 
%         [index_candidates,~]=sortrows(index_candidates,4,'ascend'); % median distance
        index_seeds=index_candidates(1:Nr_seed,:);

        % fit to obtain orientation by maximizing the completeness value,
        % constrain fitting
        tic
        Spots_pair{i}=zeros(length(index_seeds(:,1)),11);
        [indexing_final(i,:),confident_index(i),Spots_pair{i}]=fit_OR_constrain_with_local_comp_cuda(i,index_seeds,RotDet, ...
            proj_bin_bw,Spots,pos_indexing,rot_angles,S,B,Ahkl,nrhkl, ...
            thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        indexing_final0(i,:)=indexing_final(i,:);
        sprintf('Candidate %d is identified for yielding the best solution.',confident_index(i))
        if false
        if ~isempty(Spots_pair{i})
		    % fit the grain center of mass based on centers of the paired experimental spots
		    [FitOutput,fval]=fit_grain_COM(Spots_pair{i},S,B,Lsam2sou,Lsam2det,P0y,P0z,RotAxisOffset, ...
			    pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet);
		    if sqrt(sum((FitOutput(1:3)-indexing_final0(i,2:4)).^2))/VoxSize>maxD % update the voxel position, being closer to the grain COM
			    sprintf('The spot centroid difference is %.2f pixels after fitting the grain COM.',fval)
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
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
                        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
                    indexing_final(i,7)=indexing_final(i,6)/indexing_final(i,5);
                    sprintf('Update the voxel position from [%.4f %.4f %.4f] to [%.4f %.4f %.4f] mm,\n euler angles from [%.2f %.2f %.2f] to [%.2f %.2f %.2f] degrees,\n completeness from %.2f to %.2f', ...
                        pos_indexing,indexing_final(i,2:4),indexing_final0(i,9:11),indexing_final(i,9:11),indexing_final0(i,7),indexing_final(i,7))
                    pos_indexing=FitOutput(1:3);
                    pos_indices=pos2ind(i,pos_indexing,tomo_scale,VoxSize,RecVolumePixel,size(Completeness),simap_data_flag);
                    pos_list_new(i,:)=pos_indexing;
                end
            end
        end
        end
        fit_time=toc;
        sprintf('The OR calc and the fitting take %0.2f s and %0.2f s, respectively',OR_time,fit_time)
        sprintf('Indexed completeness is %.3f and Euler angles are [%.2f, %.2f, %.2f] degrees',indexing_final(i,[7 9:11]))
        
        % update
        if indexing_final(i,7)<minComp || indexing_final(i,7)<Completeness(pos_indices(1),pos_indices(2),pos_indices(3)) ...
                || indexing_final(i,8)>maxDmedian
            VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))=0.8;
            sprintf('Indexing for voxel %d failed because the completeness (%.3f) is too low (< %.3f or < %.3f) or the median D (%.2f) is too large (>%.1f pixel).', ...
                i,indexing_final(i,7),minComp,Completeness(pos_indices(1),pos_indices(2),pos_indices(3)),indexing_final(i,8),maxDmedian)
        else
            sprintf('Indexing for voxel %d is successful because the completeness (%.3f) is >= %.3f and the median D (%.2f) is <= %.1f pixels.', ...
                i,indexing_final(i,7),minComp,indexing_final(i,8),maxDmedian)
            GrainId(pos_indices(1),pos_indices(2),pos_indices(3))=i+FirstGrainID-1;
            Rodrigues(:,pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,12:14)';
            IPF001(:,pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,9:11)'./norm(indexing_final(i,9:11)); % colored by Euler angles
            EulerZXZ(:,pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,9:11)';
            VisitFlag(pos_indices(1),pos_indices(2),pos_indices(3))=1;
            Completeness(pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,7);
            Dismedian(pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,8);
            Ninter(pos_indices(1),pos_indices(2),pos_indices(3))=indexing_final(i,6);
        end
%         save(fullfile(OutputFolder,strcat(num2str(iter), 'indexing_final.mat')),'indexing_final','-v7.3');
        save(fullfile(OutputFolder,'indexing_final_temp.mat'),'indexing_final','-v7.3');
    end
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
% save(fullfile(OutputFolder,strcat(num2str(iter), 'DS_indexing.mat')),'DS_out','-v7.3');


