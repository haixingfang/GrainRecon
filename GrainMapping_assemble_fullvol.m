% postprocessing the full-volume reconstruction
% Dec 9, 2021
clear all;
close all;

SampleFlag=3;
if SampleFlag==2
    SampleName='AlCu8wt_middle_thinned_0930';
    % RecVolumePixel=[1   395;
    %                 1   399;
    %                 1   147]; % full volume
    RecVolumePixel_FOV=[63   311;
                    63   307;
                    45   115]; % effective volume masked by tomo
    OutputFolder='./AlCu_8wt_middle_thinned_0930_rec/fullvol_v3';
    fname_prefix='fullvol_v3';
elseif SampleFlag==3
    SampleName='simu_Fe';
    % RecVolumePixel=[1   300;
    %                 1   300;
    %                 1   179]; % full volume
    % RecVolumePixel_FOV=[90   210;
    %                 90   210;
    %                 1   179]; % effective volume masked by tomo, old
    RecVolumePixel_FOV=[1   160;
                1   160;
                1   240]; % effective volume masked by tomo
%    OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_v3';
%    OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_v3_laue';
%    OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_v3_simap';
%    OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_v4_laue_Icorr';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_laue_gpu_cuda_Gt_rep';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_simap_gpu_Gt';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_simap_gpu_cuda_index_compete';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_laue_old_spots_gpu_cuda_comp';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_simap_new_spots_maxD_2std_gpu_cuda_Gt';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_laue_maxD_2std_gpu_cuda_comp';
%     OutputFolder='./Fe_100um_11_11_simu_rec/fullvol_laue_new_spots_maxD_2std_gpu_cuda_comp';
%     OutputFolder='./Fe_100um_11_11_simu_rec/esrf_cluster/fullvol_simap_181projs_maxD_2std_gpu_cuda_Gt';
    OutputFolder='./Fe_100um_11_11_simu_rec/esrf_cluster/fullvol_laue_new_spots_maxD_2std_gpu_cuda_Gt';
    fname_prefix='fullvol';
elseif SampleFlag==4
    SampleName='virtual_Fe_100um_6grains';
    RecVolumePixel_FOV=[1   60;
                1   60;
                1   80]; % effective volume masked by tomo
%     OutputFolder='./virtual_Fe_100um_6grains_rec/test_round_less_scaling';
%     OutputFolder='./virtual_Fe_100um_6grains_rec/test_minus0p5_spots';
%     OutputFolder='./virtual_Fe_100um_6grains_rec/fullvol_simap_new_spots_gpu_cuda_Gt';
%     OutputFolder='./virtual_Fe_100um_6grains_rec/fullvol_simap_new_spots_gpu_cuda_Gt_esrf';
    OutputFolder='./virtual_Fe_100um_6grains_rec/test';
    fname_prefix='fullvol';
elseif SampleFlag==5
    SampleName='AlCu_8wt_middle_thinned_micro_source';
    % RecVolumePixel=[1   289;
    %                 1   313;
    %                 1   104]; % full volume
    RecVolumePixel_FOV=[28   271;
                    59   297;
                    7   87]; % effective volume masked by tomo
    OutputFolder='./AlCu_8wt_middle_thinned_micro_source_rec/fullvol_v3';
    fname_prefix='fullvol_v3';
elseif SampleFlag==6
   SampleName='AlCu_8wt_580C_1h_micro_source';
    % RecVolumePixel=[1   357;
    %                 1   343;
    %                 1   434]; % full volume
    RecVolumePixel_FOV=[95   324;
                    33   311;
                    150   430]; % effective volume masked by tomo: Z_full [16 430] but between [116 430] and [172 430] 150
    OutputFolder='./AlCu_8wt_580C_1h_micro_source_rec/fullvol_v3_linux';
    fname_prefix='fullvol_v3';
elseif SampleFlag==7
   SampleName='AlCu_8wt_580C_1h_nano_source';
    % RecVolumePixel=[1   338;
    %                 1   305;
    %                  1   425]; % full volume
    RecVolumePixel_FOV=[1   338;
                    1   305;
                    150   425]; % effective volume masked by tomo
    % OutputFolder='./AlCu_8wt_580C_1h_nano_source_rec/fullvol_no_fit';
    OutputFolder='./AlCu_8wt_580C_1h_nano_source_rec/fullvol_fitted_geo';
    fname_prefix='fullvol_v4';
else
    error('No such sample exists')
end
RecVolumePixel=RecVolumePixel_FOV;

% set up all parameters: geometry, detector, sample, reconstruction, filefolders
% SampleName='AlCu8wt_middle_thinned_0930'; % 'simu_Fe'; 'AlCu8wt_middle_thinned_0930'
setup_para;
% B=FormB(cell);
% V = cellvolume(cell); % [Angs^3]
sprintf('Tomo and spot files will be loaded from %s',FileFolder)
sprintf('Output files will be written to %s',OutputFolder)

% load tomographic volume data
% h5FileName_tomo='tomo_2022_01_20_AlCu_8wt_580C_1h_microS_vol.h5';
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
SampleVolumeDim=tomo.Dimension'.*tomo.VoxSize; % [mm]

%% discretize orientation space
% nDiv=30;
% OR=ori_division(nDiv);
% OR=ori_division_mtex('cubic'); % using mtex
OR_folder='./ori_set_hyperspherical';
OR=get_ori_set(OR_folder,sgno);

VoxSize=tomo.VoxSize(1); % [x y z] [mm]
% scale the tomo volume to have the same pixel size as RecVolume
tomo_scale=tomo;
if all(VoxSize~=tomo.VoxSize)
    tomo_scale.PhaseId=imresize3(tomo.PhaseId,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Mask=imresize3(tomo.Mask,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Dimension=size(tomo_scale.Mask);
    tomo_scale.VoxSize=[VoxSize VoxSize VoxSize];
end

% initialize
dim=size(tomo.Mask);
% load DS data
files_h5=dir([OutputFolder '/*.h5']);
files_mat=dir([OutputFolder '/*.mat']);
matFileName_prefix=[];
f_iter=[];
DS_final_flag=0;
for m=1:length(files_mat)
    matFileName=files_mat(m).name;
    if contains (matFileName, 'DS')
        is_num=~isletter(matFileName);
        if is_num(1)==1 && is_num(2)==0
            iter=str2num(matFileName(1));
        elseif is_num(1)==1 && is_num(2)==1 && is_num(3)==0
            iter=str2num(matFileName(1:2));
        elseif is_num(1)==1 && is_num(2)==1 && is_num(3)==1 && is_num(4)==0
            iter=str2num(matFileName(1:3));
        elseif contains (matFileName, 'DS_final')
            iter=0;
            DS_final_flag=1;
        else
            iter=0;
        end
        matFileName_prefix=[matFileName_prefix;iter];
        if iter>0
            load(fullfile(OutputFolder,matFileName));
            f_iter(iter)=length(find(DS_out.GrainId>0))/length(find(DS_out.Mask==1)); % indexed fraction
        end
    end
end

if DS_final_flag==0
    [~,matFileName_ind]=max(matFileName_prefix);    
    matFileName=files_mat(matFileName_ind).name;
else
    matFileName = 'DS_final.mat';
end
sprintf('Loading the file %s ...',matFileName)
load(fullfile(OutputFolder,matFileName));

if ~all(size(DS_out.PhaseId)==dim)
    DS_temp=DS_out;
    DS.PhaseId=zeros(dim);
    DS.Mask=zeros(dim);
    DS.Completeness=zeros(dim);
    DS.GrainId=zeros(dim);
    DS.Rodrigues=zeros([3 dim]);
    DS.EulerZXZ=zeros([3 dim]);
    DS.IPF001=zeros([3 dim]);
    DS.Dismedian=zeros(dim)+10; % median distance [mm]
    DS.Icorr=zeros(dim);     % intended for Icorr, but leave as no-use for the moment 
    DS.VisitFlag=zeros(dim)-1;

    clear ind ind1 ind2;
    [ind(:,1), ind(:,2), ind(:,3)]=ind2sub(size(DS_temp.Mask),find(DS_temp.Mask==1));
    [ind1(:,1), ind1(:,2), ind1(:,3)]=ind2sub(size(DS_temp.GrainId),find(DS_temp.GrainId>0));
    [ind2(:,1), ind2(:,2), ind2(:,3)]=ind2sub(size(DS_temp.Completeness),find(DS_temp.Completeness>0));
    for j=1:length(ind(:,1))
        DS.Mask(RecVolumePixel(1,1)+ind(j,1)-1,RecVolumePixel(2,1)+ind(j,2)-1, ...
            RecVolumePixel(3,1)+ind(j,3)-1)=DS_temp.Mask(ind(j,1),ind(j,2),ind(j,3));
    end
    for j=1:length(ind1(:,1))
        DS.PhaseId(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.PhaseId(ind1(j,1),ind1(j,2),ind1(j,3));
        DS.Completeness(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Completeness(ind1(j,1),ind1(j,2),ind1(j,3));
        if isfield(DS_temp,'Dismedian')
            DS.Dismedian(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Dismedian(ind1(j,1),ind1(j,2),ind1(j,3));
        end
        if isfield(DS_temp,'Icorr')
            DS.Icorr(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Icorr(ind1(j,1),ind1(j,2),ind1(j,3));
        end
        DS.GrainId(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.GrainId(ind1(j,1),ind1(j,2),ind1(j,3));
        DS.Rodrigues(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Rodrigues(:,ind1(j,1),ind1(j,2),ind1(j,3));
        DS.EulerZXZ(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.EulerZXZ(:,ind1(j,1),ind1(j,2),ind1(j,3));
        DS.IPF001(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.IPF001(:,ind1(j,1),ind1(j,2),ind1(j,3));

        DS.VisitFlag(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
            RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.VisitFlag(ind1(j,1),ind1(j,2),ind1(j,3));
    end
    % this has to put back
    RecVolumePixel=[1   size(tomo_scale.PhaseId,1);
                    1   size(tomo_scale.PhaseId,2);
                    1   size(tomo_scale.PhaseId,3)]; % full volume
else
    DS=DS_out;
end

% to remove the indexed voxels lying outside the tomo mask
ind_mask=[];
[ind_mask(:,1),ind_mask(:,2),ind_mask(:,3)]=ind2sub(size(DS.Mask),find(DS.Mask==1));
ind_GID=[];
[ind_GID(:,1),ind_GID(:,2),ind_GID(:,3)]=ind2sub(size(DS.GrainId),find(DS.GrainId>0));
ind_toremove=setdiff(ind_GID,ind_mask,'rows');
if ~isempty(ind_toremove)
    for i=1:length(ind_toremove(:,1))
        DS.Completeness(ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=0;
        DS.GrainId(ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=0;
        DS.Rodrigues(:,ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=[0 0 0]';
        DS.EulerZXZ(:,ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=[0 0 0]';
        DS.IPF001(:,ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=[0 0 0]';
        DS.Dismedian(ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=20;
        DS.VisitFlag(ind_toremove(i,1),ind_toremove(i,2),ind_toremove(i,3))=-1;
    end
end
f_indexed=length(find(DS.GrainId>0))/length(find(DS.Mask==1)); % indexed fraction

% merge regions that have smaller misorientation than pre-defined threshold value
mtex_avail=0;  % availability of mtex toolbox, 1: yes; 0: no (default).
if mtex_avail~=0
    cs = crystalSymmetry('cubic'); % 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'trigonal'
else
    cs = [];
end
min_misori = 0.5; % recommend to be 0.5 [deg]
[DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_identify_grains(DS,mtex_avail,cs,min_misori,proj_bin_bw,Spots, ...
                    rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                    RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,minComp,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                    tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
% [DS_merge]=revise_single_unindexed_voxel(DS_merge0);
sprintf('There are %d voxels with C > 1 and C_max = %.3f.',length(find(DS_merge.Completeness>1)),max(DS_merge.Completeness(:)))
DS_merge.Completeness(DS_merge.Completeness>1)=1;

% write output as h5 file and generate dream3d and xdmf files for visualization
Center=tomo.Center; % [mm]
CenterShift=[];
PhaseNo=[];
PhaseName=atomparam.name;
fname_prefix=[FileFolder(3:end) sprintf('_v4_f0p%d',round(f_indexed*100))]
ProjectName=fname_prefix;
regenerate_IPF=1;
ResultsFolder=OutputFolder;
[DS_new]=GrainMapping_writer(DS_merge,Center,CenterShift,PhaseNo,PhaseName, ...
    VoxSize,RecVolumePixel,tomo_scale.Dimension,ResultsFolder,fname_prefix,ProjectName,regenerate_IPF);

single_voxel_ID=find(DS_new.nVox==1);
if ~isempty(single_voxel_ID)
    if false
        DS_merge=revise_single_voxel(DS_new,DS_merge);
        DS_new=GrainMapping_writer(DS_merge,Center,CenterShift,PhaseNo,PhaseName, ...
            VoxSize,RecVolumePixel,tomo_scale.Dimension,ResultsFolder,fname_prefix,ProjectName,regenerate_IPF);
    end
end

% basic grain information
grains=length(DS_new.SeedID); % number of grains
grainvolume=DS_new.nVox*DS_new.VoxSize(1)*DS_new.VoxSize(2)*DS_new.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS_new.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
if false
    figure;subplot(1,2,1);hist(grainsize);subplot(1,2,2);hist(DS_new.SeedComp);
end
GrainInfo=zeros(length(DS_new.SeedID),12);
% for i=[15 24 37 38 47 50 62 73 77 85 101 102 123]
for i=1:length(DS_new.SeedID)
    if DS_new.nVox(i)>0
        pos=((DS_new.Coord(i,:)+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
        if simap_data_flag==1
            pos(1:2)=-pos(1:2);
        end
        U=euler2u(DS_new.EulerZXZ(i,1)*pi/180,DS_new.EulerZXZ(i,2)*pi/180,DS_new.EulerZXZ(i,3)*pi/180);
        [Nr_simu,Nr_intersect,dis_median,SimuSpots,HittedSpots]=index_verify_v3(U,proj_bin_bw,Spots,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        GrainInfo(i,:)=[i pos DS_new.EulerZXZ(i,:) grainsize(i) Nr_intersect Nr_simu Nr_intersect/Nr_simu dis_median];
    end
end
dlmwrite(fullfile(OutputFolder,'GrainInfo.txt'),GrainInfo,'delimiter',' ');

% below is for comparison
if SampleFlag==2
%     refFileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915';
%     h5FileName='SR_DCT_AlCu6wt_middle_thinned_trans2.h5';
    refFileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915_wolfgang';
    h5FileName='SR_DCT_AlCu6wt_middle_thinned_trans.h5';
elseif SampleFlag==3
%     refFileFolder='D:\Documents\Matlab\GrainRecon\Fe_100um_11_11_simu\input_structure';
    refFileFolder='./Fe_100um_11_11_simu/input_structure';
    h5FileName='Grain100um_400_400_600_input.h5';
elseif SampleFlag==4
%     refFileFolder='D:\Documents\Matlab\LabDCT_simap\virtual_Fe_100um_6grains';
    refFileFolder='./virtual_Fe_100um_6grains';
    h5FileName='Grain100um_6grains_150_150_200_input.h5';
elseif SampleFlag==5
    refFileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915_wolfgang';
    h5FileName='SR_DCT_AlCu6wt_middle_thinned_trans.h5';    
elseif SampleFlag==6
    refFileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu8wt_580C_1h_20220125';
%     h5FileName='SR_DCT_AlCu8wt_580C_1h_trans.h5';
    h5FileName='SR_DCT_AlCu8wt_580C_1h_trans_crop.h5';
end
Prefix_Name=h5FileName(1:end-3);
DS_ref=readLabDCT(fullfile(refFileFolder,h5FileName));

if false
    if simap_data_flag==1 && SampleFlag==4
        im_mask=[];
        for i=1:size(tomo.Mask,3)
            im_mask(:,:,i)=rot90(tomo.Mask(:,:,i),2);
        end
        DS_ref.GIDvol=uint16(double(DS_ref.GIDvol).*im_mask);
        DS_ref.CompVol=DS_ref.CompVol.*im_mask;
        DS_ref.Mask=double(DS_ref.Mask).*im_mask;
        DS_ref.RodVec3D=DS_ref.RodVec3D.*permute(repmat(im_mask,1,1,1,3),[4 1 2 3]);
        DS_ref.EulerAngle=DS_ref.EulerAngle.*permute(repmat(im_mask,1,1,1,3),[4 1 2 3]);
        DS_ref.IPF001=DS_ref.IPF001.*permute(repmat(im_mask,1,1,1,3),[4 1 2 3]);
    end
end
if SampleFlag==2
    DS_ref.Center=[0.0780   0.1371    0.2176]';
    trans_pos=[-0.0219 -0.0080 -0.0273]; % for trans2, Dec 13
    trans_pos(1:2)=-trans_pos(1:2);
    DS_ref.Center=DS_ref.Center+trans_pos';
%     center_shift=FitPos(DS_ref,DS_new,PairIndex); % [mm]
%     DS_new.Center=[0.0700 0.0600 0.1750]';
%     center_shift=[-0.0244   -0.0345    0.0173]';
%     DS_new.Center=DS_new.Center+center_shift;
end
if SampleFlag==6
    DS_new0=DS_new;
        cut_slice=[249 434];
%     cut_slice=[1 186];
    DS_merge_cut.PhaseId=DS_merge.PhaseId(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.Mask=DS_merge.Mask(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.Completeness=DS_merge.Completeness(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.GrainId=DS_merge.GrainId(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.Rodrigues=DS_merge.Rodrigues(:,:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.EulerZXZ=DS_merge.EulerZXZ(:,:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.IPF001=DS_merge.IPF001(:,:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.Dismedian=DS_merge.Dismedian(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.Icorr=DS_merge.Icorr(:,:,cut_slice(1):cut_slice(2));
    DS_merge_cut.VisitFlag=DS_merge.VisitFlag(:,:,cut_slice(1):cut_slice(2));
    
    Center=tomo_scale.Center+[0 0 tomo_scale.VoxSize(3)*(cut_slice(2)-cut_slice(1))/2]';
    fname_prefix=[FileFolder(3:end) '_cut'];
    ProjectName=fname_prefix;
    [DS_new]=GrainMapping_writer(DS_merge_cut,Center,CenterShift,PhaseNo,PhaseName, ...
        VoxSize,RecVolumePixel,tomo_scale.Dimension,ResultsFolder,fname_prefix,ProjectName,regenerate_IPF);
end
[PairIndex0,Unpaired]=DS_pair_cmp(DS_ref,DS_new);
grainvolume_ref=DS_ref.nVox*DS_ref.VoxSize(1)*DS_ref.VoxSize(2)*DS_ref.VoxSize(3)*1e9; % [um^3]
grainsize_ref=2*(3*grainvolume_ref/(4*pi)).^(1/3); % equivalent diameter of grain [um]
PairIndex=PairIndex0(PairIndex0(:,10)>=10,:); % remove the small ones
sprintf('Number of paired and uniquely paired grains in reference dataset and uniquely paired grains in rec dataset (grain size >10 voxels):')
[length(PairIndex(:,1)) length(unique(PairIndex(:,1))) length(unique(PairIndex(:,2)))]
dlmwrite(fullfile(ResultsFolder,'pairindex.txt'),PairIndex,'delimiter',' ')
if ~isempty(Unpaired)
    dlmwrite(fullfile(ResultsFolder,'Unpaired.txt'),Unpaired,'delimiter',' ')
end
[mean(PairIndex(:,10)) std(PairIndex(:,10))] % grain size
[mean(PairIndex(:,12)) std(PairIndex(:,12))] % disorientation
[mean(PairIndex(:,13)) std(PairIndex(:,13))] % grain COM [pixel]
[mean(abs(PairIndex(:,10)-PairIndex(:,9))./PairIndex(:,9)) std(abs(PairIndex(:,10)-PairIndex(:,9))./PairIndex(:,9))]

figure;
subplot(1,3,1);
hist(PairIndex(:,12));
xlabel('Disorientation (^{o})');
ylabel('Count');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',16);
axis square;
subplot(1,3,2);
hist(PairIndex(:,13));
% h1 = histogram(PairIndex(:,12));
% h1.Normalization = 'probability';
xlabel('Distance for grain COM (pixel)');
ylabel('Count');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',16);
axis square;
subplot(1,3,3);
hold all;
plot(PairIndex(:,9),PairIndex(:,10),'ro','MarkerSize',10);
plot([floor(min(min(PairIndex(:,9:10)))-10) ceil(max(max(PairIndex(:,9:10)))+10)], ...
    [floor(min(min(PairIndex(:,9:10)))-10) ceil(max(max(PairIndex(:,9:10)))+10)],'k-');
xlim([floor(min(min(PairIndex(:,9:10)))-5) ceil(max(max(PairIndex(:,9:10)))+5)]);
ylim([floor(min(min(PairIndex(:,9:10)))-5) ceil(max(max(PairIndex(:,9:10)))+5)]);
xlabel('D_{ref} (\mum)');
ylabel('D_{rec} (\mum)');
box on;
set(gca,'LineWidth',1.5);
set(gca,'FontSize',16);
axis square;
grainsize_diff=(PairIndex(:,10)-PairIndex(:,9))./PairIndex(:,9);
% print(fullfile(ResultsFolder,'cmp'),'-dtiff','-r300');
% slice_ref=slice_show(DS_ref,round(DS_ref.Dimension(3)/2),0,0,1)
% slice_rec=slice_show(DS_rec,round(DS_rec.Dimension(3)/2),0,0,1)

% visulization of overlay
SliceNumber_ref=round(DS_ref.Dimension(3)/2);
SliceNumber_rec=round(DS_new.Dimension(3)/2);
compare_XY=0;
compare_XZ=1;
compare_YZ=0;
slice_ref=slice_show(DS_ref,SliceNumber_ref,compare_XY,compare_XZ,compare_YZ)
% im=double(slice_ref);
% slice_ref_flip=joinchannels('RGB',fliplr(im(:,:,1)),fliplr(im(:,:,2)), fliplr(im(:,:,3)))

if compare_XY==1
    im=DS_ref.GIDvol(:,:,SliceNumber_ref);
    Im1=(reshape(im,[DS_ref.Dimension(1) DS_ref.Dimension(2)]));
    Im2=imrotate(Im1,90);
%     Im2=(flipud((Im2)));
elseif compare_XZ==1
    im=DS_ref.GIDvol(:,SliceNumber_ref,:);
    Im1=(reshape(im,[DS_ref.Dimension(1) DS_ref.Dimension(3)]));
    Im2=imrotate(Im1,90);
elseif compare_YZ==1
    im=DS_ref.GIDvol(SliceNumber_ref,:,:);
    Im1=(reshape(im,[DS_ref.Dimension(2) DS_ref.Dimension(3)]));
end
dipshow(Im2)

slice_rec=slice_show(DS_new,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ)
% im=double(slice_rec);
% slice_rec_flip=joinchannels('RGB',fliplr(im(:,:,1)),fliplr(im(:,:,2)), fliplr(im(:,:,3)))

DS_rec=DS_new;
[GB_overlay,DistPixel]=compare_2slices(DS_rec,DS_ref,SliceNumber_ref,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ);
% [DistIm,DistPixel1]=plot_2slices_deviation(DS_ref,DS_rec,PairIndex,SliceNumber_ref,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ);
% June 9, 2022
[DistIm,DistPixel,dev_grain,pixel_dev_grain]=cmp_3D_deviation(DS_ref,DS_rec,PairIndex, ...
                        SliceNumber_ref,compare_XY,compare_XZ,compare_YZ);
% print(fullfile(OutputFolder,'DistIm'),'-dtiff','-r400');


% for SliceNumber_rec=27:33
%     slice_rec=slice_show(DS_rec,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ)
% end
save(fullfile(ResultsFolder,sprintf('results_v4_f0p%d.mat',round(f_indexed*100))),'-v7.3');
