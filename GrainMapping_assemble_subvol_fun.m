% assemble the reconstructions of subvolumes
% Nov 4, 2021
function GrainMapping_assemble_subvol_fun()

% tomo [395,399,147]
% RecVolumePixel=[1   395;
%                 1   399;
%                 1   147]; % full volume
RecVolumePixel_FOV=[63   311;
                63   307;
                45   115]; % effective volume masked by tomo
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load cell parameters;
[cell,sgno,atomparam,space_group_IT_number]=input_al_fun();
B=FormB(cell);
V = cellvolume(cell); % [Angs^3]
simap_data_flag=1; % using simap data = 1; otherwise = 0
sample_no=2; % sample number
fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no

% load experimental parameters: source, detector and distances
[Lsam2sou,Lsam2det,dety00,detz00,tilt_x,tilt_y,tilt_z,dety0,detz0, ...
	detysize,detzsize,pixelysize,pixelzsize,BeamStopY,BeamStopZ,P0y,P0z, ...
	ExpTime]=setup_exp_fun(simap_data_flag,sample_no,fitted_geo_already);
L=Lsam2sou+Lsam2det;
RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
RotDet=RotX*RotY*RotZ;

% load tomographic volume data
FileFolder='./AlCu_8wt_middle_thinned_0930';
h5Folder_tomo=FileFolder;
h5FileName_tomo='tomo_2021_09_30_AlCu_8wt_middle_thinned.h5';
tomo=get_tomo_fromh5(fullfile(h5Folder_tomo,h5FileName_tomo),simap_data_flag);
SpotsFile='Spots_AlCu_8wt_middle_thinned_0930.mat';

% load DCT images for processing and spot segmentation, get binary images
load(fullfile(FileFolder,SpotsFile));
for i=1:length(proj_bin)
    [proj_bin_bw(:,:,i),idx] = bwdist(double(proj_bin{i}));
end
rot_angles=rot_start:rot_step:rot_end;

% generate hkl for indexing
setup_hkl;

%% discretize orientation space
nDiv=30;
OR=ori_division(nDiv);
% OR=ori_division_mtex('cubic'); % using mtex

% indexing criteria
TrustComp=0.85; % trust completeness
minComp=0.45;   % minimum completeness
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
DS.VisitFlag=zeros(dim)-1;
nodata_count=0;

FirstGrainID(1)=0;
for subvolnumber=1:n^3
    OutputFolder=['./AlCu_8wt_middle_thinned_0930_rec/subvol_' num2str(subvolnumber)];
    fname_prefix=['subvol_' num2str(subvolnumber)];
    RecVolumePixel=RecVolumePixel_job{subvolnumber};
    sprintf('Get subvolume %d ...',subvolnumber)

    files_h5=dir([OutputFolder '/*.h5']);
    files_mat=dir([OutputFolder '/*.mat']);
    if ~isempty(files_mat)
        matFileName_prefix=[];
        for m=1:length(files_mat)
            matFileName=files_mat(m).name;
            is_num=~isletter(matFileName);
            if is_num(1)==1 && is_num(2)==0
                iter=str2num(matFileName(1));
            elseif is_num(1)==1 && is_num(2)==1 && is_num(3)==0
                iter=str2num(matFileName(1:2));
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
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
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
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
        iter_merge=iter_merge+1;
        Ngrain_iter(iter_merge)=Ngrain;
        if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=1
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
[DS_new]=GrainMapping_writer(DS_merge,Center,CenterShift,PhaseNo,PhaseName, ...
    VoxSize,RecVolumePixel,tomo_scale.Dimension,[OutputFolder '_assembled'],fname_prefix,ProjectName,regenerate_IPF); 
% basic grain information
grains=length(DS_new.SeedID); % number of grains
grainvolume=DS_new.nVox*DS_new.VoxSize(1)*DS_new.VoxSize(2)*DS_new.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS_new.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
if false
    figure;subplot(1,2,1);hist(grainsize);subplot(1,2,2);hist(DS_new.SeedComp);
end
save(fullfile(OutputFolder,'Assembled.mat'),'-v7.3');

