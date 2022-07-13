% assemble the volume and identify grains, output files for visulization
% May 31, 2022
function GrainMapping_assemble_vol(DataFolder,ElementName)

    
load(fullfile(DataFolder,'para.mat'));
get_para;


% set up all parameters: geometry, detector, sample, reconstruction, filefolders
switch(ElementName)
    case 'fe'
        % input_fe;
        [cell,sgno,atomparam,space_group_IT_number]=input_fe_fun();
    case 'al'
        % input_al;
        [cell,sgno,atomparam,space_group_IT_number]=input_al_fun();
    case 'ni'
        % input_ni;
        [cell,sgno,atomparam,space_group_IT_number]=input_ni_fun();
    case 'si'
        % input_si;
        [cell,sgno,atomparam,space_group_IT_number]=input_si_fun();
    otherwise
        error('this element cannot be identified!');
end
setup_exp;
L=Lsam2sou+Lsam2det;
RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
RotDet=RotX*RotY*RotZ;

hklnumber=length(unique(Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2)); % maximum is 10, recommended be at least >= 3, 4
% iter_end=10; % maximum iterations which determines the finess of gridding for indexing seeds
sprintf('hklnumber = %.0f, C_trust = %.2f, C_min = %.2f, drop_off = %.2f, minEucDis = %.3f mm, maxD = %.1f pixel and maxDmedian = %.1f pixel', ...
    hklnumber,TrustComp,minComp,drop_off,minEucDis,maxD,maxDmedian)
if simap_data_flag==1
    S=[1 0 0;0 1 0;0 0 1];
else
    S=[1 0 0;0 -1 0;0 0 1];
end
B=FormB(cell);
V = cellvolume(cell); % [Angs^3]

sprintf('Tomo and spot files will be loaded from %s',FileFolder)
sprintf('Output files will be written to %s',OutputFolder)

% load tomographic volume data
sprintf('load tomo file: %s',tomoFile)
tomo=get_tomo_fromh5(tomoFile,1);

% load DCT images for processing and spot segmentation, get binary images
sprintf('load spots file: %s',SpotsFile)
load(SpotsFile);

for i=1:length(proj_bin)
    [proj_bin_bw(:,:,i),idx] = bwdist(double(proj_bin{i}));
end
rot_angles=rot_start:rot_step:rot_end;

% generate hkl for indexing
setup_hkl;
SampleVolumeDim=tomo.Dimension'.*tomo.VoxSize; % [mm]

% scale the tomo volume to have the same pixel size as RecVolume
tomo_scale=tomo;
if all(VoxSize~=tomo.VoxSize)
    tomo_scale.PhaseId=imresize3(tomo.PhaseId,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Mask=imresize3(tomo.Mask,tomo.VoxSize(1)/VoxSize(1),'nearest');
    tomo_scale.Dimension=size(tomo_scale.Mask);
    tomo_scale.VoxSize=[VoxSize VoxSize VoxSize];
end
tomo_scale.Center=sample.center; % [mm]
tomo_scale.VoxSize=sample.VoxSize; % [mm]

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
save(fullfile(OutputFolder,sprintf('results_v4_f0p%d.mat',round(f_indexed*100))),'-v7.3');




