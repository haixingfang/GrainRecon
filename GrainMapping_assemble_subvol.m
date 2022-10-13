% assemble the reconstructions of subvolumes
% Nov 4, 2021
clear all;
close all;


% tomo [395,399,147]
% RecVolumePixel=[1   395;
%                 1   399;
%                 1   147]; % full volume
RecVolumePixel_FOV=[63   311;
                63   307;
                45   115]; % effective volume masked by tomo
%{
% define total jobs = n^3
n=3;
vol_divide_ind=[];
for i=1:3
    vol_divide_ind(i,:)=round(linspace(RecVolumePixel_FOV(i,1),RecVolumePixel_FOV(i,2),n+1));
end
for i=1:n
    for j=1:n
        for k=1:n
            ind=sub2ind([n n n],i,j,k);
            RecVolumePixel_job{ind}(1,:)=vol_divide_ind(1,i:i+1);
            RecVolumePixel_job{ind}(2,:)=vol_divide_ind(2,j:j+1);
            RecVolumePixel_job{ind}(3,:)=vol_divide_ind(3,k:k+1);
        end
    end
end
%}

% After Dec 8, 2021
% define total jobs = n(i)^3
n=[3 3 1];
vol_divide_ind=[];
for i=1:3
    vol_divide_ind{i}=round(linspace(RecVolumePixel_FOV(i,1),RecVolumePixel_FOV(i,2),n(i)+1));
end

for i=1:n(1)
    for j=1:n(2)
        for k=1:n(3)
            ind=sub2ind([n(1) n(2) n(3)],i,j,k);
            RecVolumePixel_job{ind}(1,:)=vol_divide_ind{1}(i:i+1);
            RecVolumePixel_job{ind}(2,:)=vol_divide_ind{2}(j:j+1);
            RecVolumePixel_job{ind}(3,:)=vol_divide_ind{3}(k:k+1);
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load cell parameters;
input_al;
B=FormB(cell);
V = cellvolume(cell); % [Angs^3]
simap_data_flag=1; % using simap data = 1; otherwise = 0
sample_no=2; % sample number
fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no

% load experimental parameters: source, detector and distances
setup_exp;
L=Lsam2sou+Lsam2det;
RotDet=get_det_R(tilt_x,tilt_y,tilt_z);

% load tomographic volume data
FileFolder='./AlCu_8wt_middle_thinned_0930';
h5Folder_tomo=FileFolder;
h5FileName_tomo='tomo_2021_09_30_AlCu_8wt_middle_thinned.h5';
tomo=get_tomo_fromh5(fullfile(h5Folder_tomo,h5FileName_tomo),simap_data_flag);

% load DCT images for processing and spot segmentation, get binary images
load(fullfile(FileFolder,'Spots_AlCu_8wt_middle_thinned_0930.mat'));
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

% indexing criteria
TrustComp=0.85; % trust completeness
minComp=0.35;   % minimum completeness
drop_off=0.02;  % drop-off value for growing indexed region
minEucDis=0.5*mean([pixelysize pixelzsize])*1; % minimum tolerate euclidien distance [mm]
maxD=3; % maximum acceptable completeness weighted center difference, it needs update if larger than this value [pixel]
if simap_data_flag==1
    S=[1 0 0;0 1 0;0 0 1];
else
    S=[1 0 0;0 -1 0;0 0 1];
end

VoxSize=tomo.VoxSize(1); % [x y z] [mm]
% scale the tomo volume to have the same pixel size as RecVolume
tomo_scale=tomo;
if all(VoxSize~=tomo.VoxSize)
    tomo_scale.PhaseId=imresize3(tomo.PhaseId,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Mask=imresize3(tomo.Mask,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Dimension=size(tomo_scale.Mask);
    tomo_scale.VoxSize=[VoxSize VoxSize VoxSize];
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
dim=size(tomo.Mask);
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
nodata_count=0;

FirstGrainID(1)=0;
for subvolnumber=1:length(RecVolumePixel_job)
    OutputFolder=['./AlCu_8wt_middle_thinned_0930_rec/subvol_' num2str(subvolnumber)];
    fname_prefix=['subvol_' num2str(subvolnumber)];
    RecVolumePixel=RecVolumePixel_job{subvolnumber};
    sprintf('Get subvolume %d ...',subvolnumber)

    files_h5=dir([OutputFolder '/*.h5']);
    files_mat=dir([OutputFolder '/*.mat']);
%     if ~isempty(files_h5)
%         h5FileName=files_h5.name;
%         DS_temp=readLabDCT(fullfile(OutputFolder,h5FileName)); % read h5 file that is exported by LabDCT [X*Y*Z]
%         
%         clear ind ind1;
%         [ind(:,1), ind(:,2), ind(:,3)]=ind2sub(size(DS_temp.Mask),find(DS_temp.Mask==1));
%         [ind1(:,1), ind1(:,2), ind1(:,3)]=ind2sub(size(DS_temp.GIDvol),find(DS_temp.GIDvol>0));
%         f_indexed(subvolnumber)=length(ind1(:,1))/length(ind(:,1)); % indexed fraction
%         DS.Mask(RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.Mask;
%         DS.Completeness(RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.CompVol;
%         DS.PhaseId(DS.GrainId>0)=1;
%         DS.VisitFlag(DS.GrainId>0)=1;
%         for j=1:length(ind1(:,1))
%             DS.GrainId(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
%                 RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.GIDvol(ind1(j,1),ind1(j,2),ind1(j,3))+FirstGrainID(subvolnumber);
%         end        
%         DS.Rodrigues(:,RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.RodVec3D;
%         DS.EulerZXZ(:,RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.EulerAngle;
%         DS.IPF001(:,RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.IPF001;       
%     else
    if ~isempty(files_mat)
        matFileName_prefix=[];
        for m=1:length(files_mat)
            matFileName=files_mat(m).name;
            is_num=~isletter(matFileName);
            if is_num(1)==1 && is_num(2)==0
                iter=str2num(matFileName(1));
            elseif is_num(1)==1 && is_num(2)==1 && is_num(3)==0
                iter=str2num(matFileName(1:2));
            elseif is_num(1)==1 && is_num(2)==1 && is_num(3)==1 && is_num(4)==0
                iter=str2num(matFileName(1:3));
            end
            matFileName_prefix=[matFileName_prefix;iter];
        end
        [~,matFileName_ind]=max(matFileName_prefix);    
        matFileName=files_mat(matFileName_ind).name;
        sprintf('Loading the file %s ...',matFileName)
        load(fullfile(OutputFolder,matFileName));
        DS_temp=DS_out;
        
        clear ind ind1 ind2;
        [ind(:,1), ind(:,2), ind(:,3)]=ind2sub(size(DS_temp.Mask),find(DS_temp.Mask==1));
        [ind1(:,1), ind1(:,2), ind1(:,3)]=ind2sub(size(DS_temp.GrainId),find(DS_temp.GrainId>0));
        [ind2(:,1), ind2(:,2), ind2(:,3)]=ind2sub(size(DS_temp.Completeness),find(DS_temp.Completeness>0));
        why_id=[];
        for m=1:length(ind1(:,1))
            rod=DS_temp.Rodrigues(:,ind1(m,1),ind1(m,2),ind1(m,3))';
            if all(rod==0)
                why_id=[why_id;ind1(m,:)];
            end
        end
        f_indexed(subvolnumber)=length(ind1(:,1))/length(ind(:,1)); % indexed fraction
%         DS.PhaseId(RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.PhaseId;
%         DS.Mask(RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.Mask;
%         DS.Completeness(RecVolumePixel(1,1):RecVolumePixel(1,2), ...
%             RecVolumePixel(2,1):RecVolumePixel(2,2),RecVolumePixel(3,1):RecVolumePixel(3,2))=DS_temp.Completeness;
        for j=1:length(ind1(:,1))
            DS.PhaseId(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.PhaseId(ind1(j,1),ind1(j,2),ind1(j,3));
            DS.Mask(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Mask(ind1(j,1),ind1(j,2),ind1(j,3));
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
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.GrainId(ind1(j,1),ind1(j,2),ind1(j,3))+FirstGrainID(subvolnumber);
            DS.Rodrigues(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.Rodrigues(:,ind1(j,1),ind1(j,2),ind1(j,3));
            DS.EulerZXZ(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.EulerZXZ(:,ind1(j,1),ind1(j,2),ind1(j,3));
            DS.IPF001(:,RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.IPF001(:,ind1(j,1),ind1(j,2),ind1(j,3));
            
            DS.VisitFlag(RecVolumePixel(1,1)+ind1(j,1)-1,RecVolumePixel(2,1)+ind1(j,2)-1, ...
                RecVolumePixel(3,1)+ind1(j,3)-1)=DS_temp.VisitFlag(ind1(j,1),ind1(j,2),ind1(j,3));
        end
        clear ind1;
        [ind1(:,1), ind1(:,2), ind1(:,3)]=ind2sub(size(DS.GrainId),find(DS.GrainId>0));
        why_id=[];
        for m=1:length(ind1(:,1))
            rod=DS.Rodrigues(:,ind1(m,1),ind1(m,2),ind1(m,3))';
            if all(rod==0)
                why_id=[why_id;ind1(m,:)];
                error('wrong assignment');
            end
        end
    else
        sprintf('No data is found for subvolume %d ...',subvolnumber)
        nodata_count=nodata_count+1;
    end
    FirstGrainID(subvolnumber+1)=max(DS.GrainId(:));
end

% this has to put back
RecVolumePixel=[1   size(tomo_scale.PhaseId,1);
                1   size(tomo_scale.PhaseId,2);
                1   size(tomo_scale.PhaseId,3)]; % full volume

% merge regions that have smaller misorientation than pre-defined threshold value
mtex_avail=0;  % availability of mtex toolbox, 1: yes; 0: no (default).
if mtex_avail~=0
    cs = crystalSymmetry('cubic'); % 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'trigonal'
else
    cs = [];
end
min_misori = 0.5; % recommend to be 0.5 [deg]
[DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness(DS,mtex_avail, ...
    cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
iter_merge=1;
Ngrain_iter(iter_merge)=Ngrain;
iter_flag=1;  % iterative merging
if iter_flag==1
    stop_iter=0;
    while stop_iter~=1
        [DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness(DS_merge,mtex_avail, ...
            cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
        iter_merge=iter_merge+1;
        Ngrain_iter(iter_merge)=Ngrain;
        if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=0
            stop_iter=1;
        end
    end
end

% write output as h5 file and generate dream3d and xdmf files for visualization
Center=tomo.Center; % [mm]
CenterShift=[];
PhaseNo=[];
PhaseName=atomparam.name;
fname_prefix=FileFolder(3:end);
ProjectName=fname_prefix;
regenerate_IPF=1;
ResultsFolder=[OutputFolder '_assemble_drop0p05_OR2deg_angle_1_minComp_0p3_v3'];
[DS_new]=GrainMapping_writer(DS_merge,Center,CenterShift,PhaseNo,PhaseName, ...
    VoxSize,RecVolumePixel,tomo_scale.Dimension,ResultsFolder,fname_prefix,ProjectName,regenerate_IPF); 
% basic grain information
grains=length(DS_new.SeedID); % number of grains
grainvolume=DS_new.nVox*DS_new.VoxSize(1)*DS_new.VoxSize(2)*DS_new.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS_new.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
if false
    figure;subplot(1,2,1);hist(grainsize);subplot(1,2,2);hist(DS_new.SeedComp);
end



% below is for comparison
FileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915';
h5FileName='SR_DCT_AlCu6wt_middle_thinned_trans2.h5';
Prefix_Name=h5FileName(1:end-3);
DS_ref=readLabDCT(fullfile(FileFolder,h5FileName));
DS_ref.Center=[0.0780   0.1371    0.2176]';
trans_pos=[-0.0219 -0.0080 -0.0273]; % for trans2, Dec 13
trans_pos(1:2)=-trans_pos(1:2);
DS_ref.Center=DS_ref.Center+trans_pos';
[PairIndex0,Unpaired]=DS_pair_cmp(DS_ref,DS_new);
grainvolume_ref=DS_ref.nVox*DS_ref.VoxSize(1)*DS_ref.VoxSize(2)*DS_ref.VoxSize(3)*1e9; % [um^3]
grainsize_ref=2*(3*grainvolume_ref/(4*pi)).^(1/3); % equivalent diameter of grain [um]
PairIndex=PairIndex0(PairIndex0(:,10)>=10,:); % remove the small ones
sprintf('Number of paired and uniquely paired grains in reference dataset and uniquely paired grains in rec dataset:')
[length(PairIndex(:,1)) length(unique(PairIndex(:,1))) length(unique(PairIndex(:,2)))]
dlmwrite(fullfile(ResultsFolder,'pairindex.txt'),PairIndex,'delimiter',' ')
if ~isempty(Unpaired)
    dlmwrite(fullfile(ResultsFolder,'Unpaired.txt'),Unpaired,'delimiter',' ')
end
[mean(PairIndex(:,10)) std(PairIndex(:,10))] % grain size
[mean(PairIndex(:,12)) std(PairIndex(:,12))] % disorientation
[mean(PairIndex(:,13)) std(PairIndex(:,13))] % grain COM [pixel]
[mean((PairIndex(:,10)-PairIndex(:,9))./PairIndex(:,9)) std((PairIndex(:,10)-PairIndex(:,9))./PairIndex(:,9))]
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
% print('cmp','-dtiff','-r300');
% slice_ref=slice_show(DS_ref,round(DS_ref.Dimension(3)/2),0,0,1)
% slice_rec=slice_show(DS_rec,round(DS_rec.Dimension(3)/2),0,0,1)

% visulization of overlay
DS_rec=DS_new;
SliceNumber_ref=100;
SliceNumber_rec=100;
compare_XY=0;
compare_XZ=1;
compare_YZ=0;
slice_ref=slice_show(DS_ref,SliceNumber_ref,compare_XY,compare_XZ,compare_YZ)
slice_rec=slice_show(DS_rec,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ)
[GB_overlay,DistPixel1]=compare_2slices(DS_ref,DS_rec,SliceNumber_ref,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ);
for SliceNumber_rec=27:33
    slice_rec=slice_show(DS_rec,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ)
end

save(fullfile(ResultsFolder,'results.mat'),'-v7.3');

