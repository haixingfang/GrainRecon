% calculate the completeness for all voxels in the DS file
% gpu calculation
function DS_out = Forward_calc_comp_allvoxel(DS,DataFolder,ElementName,trans_pos)

if nargin<4
    trans_pos=[0 0 0];
end
load(fullfile(DataFolder,'para.mat'));
get_para;
% set up all parameters: geometry, detector, sample, reconstruction, filefolders
switch(ElementName)
    case 'fe'
        input_fe;
        [cell,sgno,atomparam,space_group_IT_number]=input_fe_fun();
    case 'al'
        input_al;
        [cell,sgno,atomparam,space_group_IT_number]=input_al_fun();
    case 'ni'
        input_ni;
        [cell,sgno,atomparam,space_group_IT_number]=input_ni_fun();
    case 'si'
        input_si;
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
if simap_data_flag==1
    S=[1 0 0;0 1 0;0 0 1];
else
    S=[1 0 0;0 -1 0;0 0 1];
end

hklnumber=length(unique(Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2)); % maximum is 10, recommended be at least >= 3, 4
% iter_end=10; % maximum iterations which determines the finess of gridding for indexing seeds
sprintf('hklnumber = %.0f, C_trust = %.2f, C_min = %.2f, drop_off = %.2f, minEucDis = %.3f mm, maxD = %.1f pixel and maxDmedian = %.1f pixel', ...
    hklnumber,TrustComp,minComp,drop_off,minEucDis,maxD,maxDmedian)
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

% set to be consistent with the size of the whole volume
RecVolumePixel(:,1)=[1 1 1]';
RecVolumePixel(:,2)=size(DS.GIDvol)';
for i=1:3
    tomo.Center(i)=tomo.Center(i)-trans_pos(i);
end

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

% basic grain information
grains=length(DS.SeedID); % number of grains
grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
Rsample=sqrt(SampleVolumeDim(1)^2+SampleVolumeDim(2)^2)/2; % [mm]

if gpuDeviceCount( "available" )
    DS_out = Forward_calc_comp_allvoxel_fun(DS,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00, ...
                            P0y,P0z,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                            RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxD,maxDmedian);
else
    warning("GPU is required for the computation! NO processing has been done!")
    DS_out=DS;
end

