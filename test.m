% GrainMapping for lab-base diffraction contrast tomography
% July 8, 2021
% Copyright by Haixing Fang, haixing.fang@grenoble-inp.fr; haixingfang868@gmail.com
% indexing grain orientations based on forward simulation
% shape reconstruction by growing neighboring voxels around indexed seeding voxels
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% steps for grain mapping
% step 1: get geometry information using get_geometry.m and input those information to setup_exp.m
% step 2: get spots info
%         if spot segmented images have been obtained, get spot info from binarized images using get_spots.m and set load_spot = 1;
%         if DCT projections are not segmented yet, continue this program by setting load_spot = 0
% step 3: get absorption tomo, run get_tomo_slices_create_h5.m
% step 4: check other parameters, e.g. RecBox, RecCenter, RecVoxelSize, S etc;
%         note the coordinate differences betweeen DTU and Simap datasets
% step 5: run GrainMapping.m or GrainMapping_simap.m depending where the dataset is obtained
% step 5-1: input crytallographic parameters, experimental geometry, tomo h5 and spots info
% step 5-2: setup_hkl
% step 5-3: define reconstruction info: voxel size, volume, center and reconstruction parameters
% step 5-4: indexing and growing process (decoupled)
% step 5-5: merging process to define grains
% step 5-6: write out h5 ouput for DS
% step 5-7: fitting geometry for redoing the whole reconstruction
% step 5-8: end, when there is no need to redo the reconstruction given that good fitting parameters have been implemented
% version 2.0, Sep 24, 2021
% version 3.0, Dec 3, 2021
% version 4.0, Feb 10, 2022 (remove assigned spots from the spot list)
% version 5.0, April 8, 2022, compute by cpu, gpu_cuda_Gt, gpu_cuda_comp
% function GrainMapping_fun(OutputFolder,fname_prefix,RecVolumePixel,SampleName,compute_opt)
%for testing
clear all;
% SampleName='virtual_Fe_100um_6grains';
% RecVolumePixel=[1   60;
%             1   60;
%             1   80];
% OutputFolder='./virtual_Fe_100um_6grains_rec/test';
% fname_prefix='fullvol';

% SampleName='AlCu_8wt_580C_1h_nano_source';
% % RecVolumePixel=[1   338;
% %                 1   305;
% %                  1   425]; % full volume
% RecVolumePixel=[1   338;
%                 1   305;
%                 150   425]; % effective volume masked by tomo
% OutputFolder='./AlCu_8wt_580C_1h_nano_source_rec/fullvol_fitted_geo_gpu_test';
% fname_prefix='fullvol_v4';

SampleName='simu_Fe';
RecVolumePixel=[1   160;
                1   160;
                1   240];
OutputFolder='./Examples/Fe_100um_11_11_simu_rec/test';
fname_prefix='fullvol';

computation_options = [{'cpu'}, {'gpu_cuda_Gt'}, {'gpu_cuda_comp'}];
compute_opt = computation_options{3};

sprintf('Computational choice: %s',compute_opt)

% set up all parameters: geometry, detector, sample, reconstruction, filefolders
% SampleName='AlCu8wt_middle_thinned_0930'; % 'simu_Fe'; 'AlCu8wt_middle_thinned_0930'
setup_para;
sprintf('Tomo and spot files will be loaded from %s',FileFolder)
sprintf('Output files will be written to %s',OutputFolder)

% load tomographic volume data
sprintf('load tomo file: %s',fullfile(h5Folder_tomo,h5FileName_tomo))
tomo=get_tomo_fromh5(fullfile(h5Folder_tomo,h5FileName_tomo),1);

% load DCT images for processing and spot segmentation, get binary images
sprintf('load spots file: %s',fullfile(FileFolder,SpotsFile))
load(fullfile(FileFolder,SpotsFile));

for i=1:length(proj_bin)
    [proj_bin_bw(:,:,i),idx] = bwdist(double(proj_bin{i}));
end
rot_angles=rot_start:rot_step:rot_end;

% generate hkl for indexing
setup_hkl;

sprintf('Geometry: Lss = %.2f mm, Lsd = %.2f mm, dety00 = %.2f mm, detz00 = %.2f mm,\n tilt_xyz = [%.2f %.2f %.2f] degrees', ...
    Lsam2sou, Lsam2det, dety00, detz00, tilt_x, tilt_y, tilt_z)
rng('shuffle');
%% Define the reconstruction volume
% figure;view0=orthosliceViewer(tomo.PhaseId,'Colormap',parula(256));
SampleVolumeDim=tomo.Dimension'.*tomo.VoxSize; % [mm]
% check orthoviewer to get the ROI
% note: in orthoview [x y] is reversed from conventionally defined system,
% should be switched
tomo_FOV_range=[1 tomo.Dimension(1);1 tomo.Dimension(2);1 tomo.Dimension(3)]; % [pixel]
VoxSize=tomo.VoxSize(1); % [x y z] [mm]
tomo_FOV_center=(tomo_FOV_range(:,2)-tomo_FOV_range(:,1))/2+tomo_FOV_range(:,1)-0.5; % [pixel]
RecVol_origin=(-tomo.Dimension'/2+tomo_FOV_center).*tomo.VoxSize; % [mm]
RecVol_origin(1:2)=-RecVol_origin(1:2); % add on Oct 7, 2021
RecVol_shift=RecVol_origin./VoxSize; % shift with respect to the tomo volume [pixel]
for i=1:3
    RecVolume(i,1)=(RecVolumePixel(i,1)-1)*VoxSize-SampleVolumeDim(i)/2;
    RecVolume(i,2)=RecVolumePixel(i,2)*VoxSize-SampleVolumeDim(i)/2;  
    RecBox(i)=-(RecVolume(i,1)-RecVol_origin(i))*2; % [mm]
end
for i=1:3
    if RecVolumePixel(i,1)==0
        RecVolumePixel(i,1)=1;
        RecVolumePixel(i,2)=RecVolumePixel(i,2)+1;
    end
end
dim=RecVolumePixel(:,2)-RecVolumePixel(:,1)+1;
dim=dim';
sprintf('The reconstructed volume has dimensions of %d*%d*%d pixel',dim)

if dim(1)>tomo.Dimension(1) || dim(2)>tomo.Dimension(2) || dim(3)>tomo.Dimension(3) 
    error('Error: reconstructed volume is larger than sample volume, please set a smaller volume for reconstruction.')
elseif dim(3)==1
    error('Error: please set the volume_z >= 2 * VoxSize.')
end

% scale the tomo volume to have the same pixel size as RecVolume
tomo_scale=tomo;
if all(VoxSize~=tomo.VoxSize)
    tomo_scale.PhaseId=imresize3(tomo.PhaseId,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Mask=imresize3(tomo.Mask,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Dimension=size(tomo_scale.Mask);
    tomo_scale.VoxSize=[VoxSize VoxSize VoxSize];
end
% figure;view1=orthosliceViewer(tomo_scale.PhaseId,'Colormap',parula(256));
if RecVolumePixel(1,2)>tomo_scale.Dimension(1) || RecVolumePixel(2,2)>tomo_scale.Dimension(2) ...
        || RecVolumePixel(3,2)>tomo_scale.Dimension(3) || ~all(RecVolumePixel(:,1)>0)
    error('Error: reconstructed volume is not properly set, please reset it.')
end

% initialize information for the reconstructed volume
DS.PhaseId=zeros(dim);
DS.Mask=zeros(dim);
DS.Completeness=zeros(dim);
DS.GrainId=zeros(dim);
DS.Rodrigues=zeros([3 dim]);
DS.EulerZXZ=zeros([3 dim]);
DS.IPF001=zeros([3 dim]);
DS.Dismedian=zeros(dim)+10; % median distance [mm]
DS.Ninter=zeros(dim);    % number of intersected spots
DS.Icorr=zeros(dim);     % intended for Icorr, but leave as no-use for the moment 
DS.VisitFlag=zeros(dim)-1; % flag for visiting indexing, -1: non-visible; 0: not visited; 1: visited
clear allVoxel_indices;
[allVoxel_indices(:,1), allVoxel_indices(:,2), allVoxel_indices(:,3)]=ind2sub(dim, ...
    find(tomo_scale.Mask(RecVolumePixel(1,1):RecVolumePixel(1,2),RecVolumePixel(2,1):RecVolumePixel(2,2), ...
    RecVolumePixel(3,1):RecVolumePixel(3,2))==1));
for i=1:length(allVoxel_indices(:,1))
    DS.Mask(allVoxel_indices(i,1),allVoxel_indices(i,2),allVoxel_indices(i,3))=1;
    DS.PhaseId(allVoxel_indices(i,1),allVoxel_indices(i,2),allVoxel_indices(i,3))=1;
    DS.VisitFlag(allVoxel_indices(i,1),allVoxel_indices(i,2),allVoxel_indices(i,3))=0;
    DS.Dismedian(allVoxel_indices(i,1),allVoxel_indices(i,2),allVoxel_indices(i,3))=20;
end
allVoxel_indices=allVoxel_indices+RecVolumePixel(:,1)'-1;

if isempty(allVoxel_indices)
    error('Error: please select a sample region containing grains to be reconstructed.');
end
allVoxel_indices=allVoxel_indices(randperm(length(allVoxel_indices(:,1))),:);
total_voxel=length(allVoxel_indices(:,1));

% generate seeding positions
pos_seed=generate_seeding_indexing_pos(1,allVoxel_indices,tomo_scale,VoxSize,simap_data_flag);

%% discretize orientation space
% nDiv=30;
% OR=ori_division(nDiv);
% OR=ori_division_mtex('cubic'); % using mtex
OR_folder='./ori_set_hyperspherical';
if strcmp(compute_opt, 'gpu_cuda_comp')
    OR=get_ori_set(OR_folder,sgno,'1')
else
    OR=get_ori_set(OR_folder,sgno)
end
vis_flag=0;
if vis_flag==1
    cs = crystalSymmetry('cubic');
    OR_ori=orientation(rotation(quaternion(OR.q(:,1),OR.q(:,2),OR.q(:,3),OR.q(:,4))),cs);
    figure;
    plot(OR_ori,'axisAngle','all');
end

%%%%%%%%%%%%%%%%%% grain mapping engine
SpotsForIndexing_select=0; % select spots which above certain size for indexing
if SpotsForIndexing_select==1
    for i=1:length(Spots)
        SpotsForIndex{i}=Spots{i}(Spots{i}(:,11)>20,:); % e.g. sizes should > 100 pixels
    end
else
    SpotsForIndex=Spots;
end

total_run_time=0; % total running time [s]
iter=0;
stop_grain_mapping=0;
indexed_voxel_fraction_min=0.999; % minimum allowed fraction of indexed voxels
% indexed_voxel_fraction_min=20000/total_voxel; % minimum allowed fraction of indexed voxels
FirstGrainID(1)=1;

indexing_final_all=[];
pos_seed_all=[];
Nr_indexed_voxel_all=[];
Nr_seed=50;
% save parameters
save_para(Lsam2sou,Lsam2det,P0y,P0z,dety00,detz00,tilt_x,tilt_y,tilt_z,RotDet, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
    drop_off,maxD,TrustComp,minComp,minEucDis,Ahkl,OR,VoxSize,RecVolumePixel, ...
    FileFolder,OutputFolder,simap_data_flag,sample_no,fitted_geo_already,rot_angles, ...
    tomo_scale,h5Folder_tomo,h5FileName_tomo,SpotsFile,maxDmedian);
% save(fullfile(OutputFolder,'recon_init.mat'),'-v7.3');

if strcmp(compute_opt, 'cpu')
    if isempty(gcp('nocreate'))
        parpool('local',maxNumCompThreads);
    end
elseif strcmp(compute_opt, 'gpu_cuda_Gt')
    fprintf(['Note: you may need to modify the pre-assigned size for "float Gt_match" in "cuda_forward_Gt.cu".\n' ...
    'It is recomended to set it as "float Gt_match[nr_spots*7]" where nr_spots = nr_projs*maximum spots per projection. \n']);
    nr_spots = 0;
    for i=1:length(Spots)
        nr_spots = nr_spots + length(Spots{i}(:,1));
    end
    sprintf('There are in total %d spots.', nr_spots)
end

if ~exist('iter_end','var')
    iter_end = 10; % define the gridding level
elseif iter_end > 10
    iter_end=10;
end
check_fit_all=0;
if check_fit_all==1
    sprintf("Grain mapping engine will check the need of fit_all mode for seeds.")
end
while stop_grain_mapping~=1
    iter=iter+1;
    % generate seeding positions for indexing, April 14, 2022
%     pos_seed = generate_uniform_seeding_pos(iter,RecVolumePixel,tomo_scale,VoxSize,simap_data_flag);

    tStart=tic;
    sprintf('Grain mapping engine starts with an iterative attempt of %d ...',iter)
    sprintf('There are %d / %d unvisited voxels and %d / %d = %.3f voxels to be indexed ...', ...
        length(find(DS.VisitFlag==0)),total_voxel,total_voxel-length(find(DS.GrainId>0)), ...
        total_voxel,(total_voxel-length(find(DS.GrainId>0)))/total_voxel)
    switch (compute_opt)
        case 'cpu'
             [indexing_final,pos_seed_new,DS_out,Nr_indexed_voxel]=GrainMapping_engine_cpu(pos_seed,OR,RotDet, ...
                SpotsForIndex,proj_bin_bw,rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                hkl_family,hkl_family_square,d_possible,Glen_possible,thetamax,lambda_min,lambda_max, ...
                Lsam2sou,Lsam2det,TrustComp,minComp,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                drop_off,maxD,DS,tomo_scale,VoxSize,RecVolumePixel,FirstGrainID(iter),simap_data_flag,maxDmedian,Nr_seed,OutputFolder,iter);
        case 'gpu_cuda_Gt'
             [indexing_final,pos_seed_new,DS_out,Nr_indexed_voxel]=GrainMapping_engine_gpu_Gt(pos_seed,OR,RotDet, ...
                SpotsForIndex,proj_bin_bw,rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                hkl_family_square,d_possible,thetamax,lambda_min,lambda_max, ...
                Lsam2sou,Lsam2det,TrustComp,minComp,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                drop_off,maxD,DS,tomo_scale,VoxSize,RecVolumePixel,FirstGrainID(iter),simap_data_flag, ...
                maxDmedian,Nr_seed,OutputFolder,iter,check_fit_all);
         case 'gpu_cuda_comp'
             [indexing_final,pos_seed_new,DS_out,Nr_indexed_voxel]=GrainMapping_engine_gpu_comp(pos_seed,OR,RotDet, ...
                SpotsForIndex,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,thetamax,lambda_min,lambda_max, ...
                Lsam2sou,Lsam2det,TrustComp,minComp,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                drop_off,maxD,DS,tomo_scale,VoxSize,RecVolumePixel,FirstGrainID(iter), ...
                simap_data_flag,maxDmedian,Nr_seed,OutputFolder,iter,check_fit_all);
    end
    iter_time=toc(tStart);
    sprintf('To complete iteration %d takes %.2f s and successful indexing rate = %.3f.',iter,iter_time, ...
        length(find(indexing_final(:,7)>=minComp))/length(find(indexing_final(:,7)>0)))
    save(fullfile(OutputFolder,[strcat('pos',num2str(iter)) '.mat']),'pos_seed','pos_seed_new');
    %%%%%%%%%%%%%%%% store intermediate data
    if length(find(indexing_final(:,7)>=minComp))/length(find(indexing_final(:,7)>0))<1/4
        Nr_seed=500;
        sprintf("Nr_seed switches its value to %d.", Nr_seed)
        not_remove_spot=1;
    end
    total_run_time=total_run_time+iter_time;
    indexing_final_all=[indexing_final_all;indexing_final];
    pos_seed_all=[pos_seed_all;pos_seed_new];
    Nr_indexed_voxel_all=[Nr_indexed_voxel_all Nr_indexed_voxel];   
    DS_record{iter}=DS_out;
    
    indexed_voxel_fraction(iter)=length(find(DS_out.GrainId>0))/total_voxel;
    if (iter>=2 && indexed_voxel_fraction(iter)>indexed_voxel_fraction_min && ...
               abs(indexed_voxel_fraction(iter)-indexed_voxel_fraction(iter-1))*total_voxel<=2) || ...
               (iter==1 && indexed_voxel_fraction(iter)>0.999) || (iter>=iter_end ...
               && (indexed_voxel_fraction(iter)>0.995 ||  (total_voxel-length(find(DS_out.GrainId>0)))<30000))
        stop_grain_mapping=1;
        sprintf('There are %d / %d voxels already indexed ...',length(find(DS_out.GrainId>0)),total_voxel)
    else  % continue to index unvisited voxel
        FirstGrainID(iter+1)=FirstGrainID(iter)+length(pos_seed(:,1));
%         clear un_visited_voxels;
%         [un_visited_voxels(:,1), un_visited_voxels(:,2), un_visited_voxels(:,3)]=ind2sub(size(DS_out.VisitFlag), ...
%            find(DS_out.VisitFlag>=0 & DS_out.VisitFlag<0.5));
%         [un_visited_voxels(:,1), un_visited_voxels(:,2), un_visited_voxels(:,3)]=ind2sub(size(DS_out.VisitFlag), ...
%             find(DS_out.VisitFlag==0));
        
        % update on April 24, 2022
        clear un_visited_voxels;
        if indexed_voxel_fraction(iter)>0.9
            completeness_bin1=DS_out.Completeness>=TrustComp;
            completeness_bin2=DS_out.VisitFlag==1;
            completeness_bin=completeness_bin1+completeness_bin2>0;
            completeness_bin_bw=bwdist(double(completeness_bin)); % the seeding candidates should far away from the indexed ones        
            [un_visited_voxels(:,1), un_visited_voxels(:,2), un_visited_voxels(:,3)]=ind2sub(size(DS_out.VisitFlag), ...
                find(DS_out.VisitFlag>=0 & completeness_bin_bw>=3));
        else
           [un_visited_voxels(:,1), un_visited_voxels(:,2), un_visited_voxels(:,3)]=ind2sub(size(DS_out.VisitFlag), ...
                        find(DS_out.VisitFlag>=0 & DS_out.VisitFlag<0.5));
        end
        if (~isempty(un_visited_voxels) && length(un_visited_voxels(:,1))<total_voxel*0.05) || isempty(un_visited_voxels)
	        clear un_visited_voxels;
            [un_visited_voxels(:,1), un_visited_voxels(:,2), un_visited_voxels(:,3)]=ind2sub(size(DS_out.VisitFlag), ...
                find(DS_out.VisitFlag>=0 & DS_out.VisitFlag<0.5));
        end        
        if ~isempty(un_visited_voxels)
            un_visited_voxels=un_visited_voxels+RecVolumePixel(:,1)'-1;
            pos_seed=generate_seeding_indexing_pos(iter+1,un_visited_voxels,tomo_scale,VoxSize,simap_data_flag);
        else
            stop_grain_mapping=1;
            sprintf('There are %d / %d voxels already indexed ...',length(find(DS_out.GrainId>0)),total_voxel)
        end
    end
    DS=DS_out;
    save(fullfile(OutputFolder,strcat(num2str(iter), 'DS.mat')),'DS_out','-v7.3');

%    % remove the assigned spots from the spot list
   if strcmp(compute_opt,'gpu_cuda_Gt')
       for k=1:length(SpotsForIndex)
           if isempty(SpotsForIndex)
               not_remove_spot=1;
           end
       end
       if iter>2 && not_remove_spot~=1
           if (indexed_voxel_fraction(iter-1) == 0 && indexed_voxel_fraction(iter)>0.3) ...
               || (indexed_voxel_fraction(iter-1)>0 && abs(indexed_voxel_fraction(iter)-indexed_voxel_fraction(iter-1))<0.1)
               [SpotsForIndex,~]=merge_regions_and_remove_spots(DS,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber, ...
                           hkl_square,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                           pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ,tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
    %            GrainMapping_writer(DS_merge,tomo_scale.Center,[],[],atomparam.name, ...
    %               VoxSize,RecVolumePixel,tomo_scale.Dimension,OutputFolder,[SampleName '_temp'],[SampleName '_temp'],1);
         end
       end
   end
end

