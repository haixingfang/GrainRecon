% Fit the geometry using the reconstructed map
% grain map is saved as a .h5 file
% thresholds of grains for fitting: minimum completeness (comp_min) and minimum grain size (grainsize_min)
% h5FileName contains '.h5'
function ParaFit = GrainMapping_fit_geo(DataFolder,ElementName,comp_min,grainsize_min,h5FileName)

load(fullfile(DataFolder,'para.mat'));
get_para;
% set up all parameters: geometry, detector, sample, reconstruction, filefolders
switch(ElementName)
    case 'fe'
%         input_fe;
        [cell,sgno,atomparam,space_group_IT_number]=input_fe_fun();
    case 'al'
%         input_al;
        [cell,sgno,atomparam,space_group_IT_number]=input_al_fun();
    case 'ni'
%         input_ni;
        [cell,sgno,atomparam,space_group_IT_number]=input_ni_fun();
    case 'si'
%         input_si;
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
if nargin<5
    files_h5=dir([OutputFolder '/*.h5']);
    h5FileName=files_h5(1).name;
end
Prefix_Name=h5FileName(1:end-3);
DS_new=readLabDCT(fullfile(OutputFolder,h5FileName));
% basic grain information
grains=length(DS_new.SeedID); % number of grains
grainvolume=DS_new.nVox*DS_new.VoxSize(1)*DS_new.VoxSize(2)*DS_new.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS_new.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])

% set to be consistent with the size of the whole volume
RecVolumePixel(:,1)=[1 1 1]';
RecVolumePixel(:,2)=size(DS.GIDvol)';

%%%%%% fit geometry, select grains for fitting
% select grains for fitting
if nargin<3
    comp_min=0.55;
    grainsize_min=45; % [um]
elseif nargin<4
    grainsize_min=45; % [um]
end
hittedSpots_pair=[];
hittedSpots_pair_all=[];
select_ind=find(grainsize>grainsize_min & DS_new.SeedComp>=comp_min);
if isempty(select_ind)
    [~,select_ind]=max(grainsize);
end
sprintf('%d grains of total %d grains are selected for fitting',length(select_ind),length(grainsize))
    
% select grains for fitting
DS_fit=DS_new;
eff_index=[];
for i=1:length(select_ind)
    eff_index=[eff_index;find(DS_new.SeedID==select_ind(i))];
end
DS_fit.SeedID=DS_new.SeedID(eff_index);
DS_fit.SeedComp=DS_new.SeedComp(eff_index);
DS_fit.RodVec=DS_new.RodVec(eff_index,:);
DS_fit.nVox=DS_new.nVox(eff_index);
DS_fit.Coord=DS_new.Coord(eff_index,:);
DS_fit.EulerZXZ=DS_new.EulerZXZ(eff_index,:);

tic
Rsample=sqrt(SampleVolumeDim(1)^2+SampleVolumeDim(2)^2)/2; % [mm]
if length(DS_new.SeedID)>length(DS_fit.SeedID)
%     rot_angles_simu=[rot_start:3*rot_step:rot_end-180];
    rot_angles_simu=0;
    [SpotNr_simu_all,SpotNr_obs_all,SpotsPair_all,GrainIndex_all,SubGrain_all,~,SmallGrID]=Forward_simu_spots_exp(DS_new, ...
        Rsample,RecVolumePixel,tomo_scale,ExpTime,atomparam,proj,Spots,rot_start,rot_step, ...
        S,B,Ahkl,nrhkl,hkl_square,Energy,lambda,V,K1,I0E,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det, ...
        dety00,detz00,P0y,P0z,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        simap_data_flag,strcat(OutputFolder,'/grains_all'),[rot_start:3*rot_step:rot_end-180]);
end
    rot_angles_simu=[rot_start:2*rot_step:rot_end];
    [SpotNr_simu,SpotNr_obs,SpotsPair,GrainIndex,SubGrain,rot_angles_calc,~]=Forward_simu_spots_exp(DS_fit, ...
        Rsample,RecVolumePixel,tomo_scale,ExpTime,atomparam,proj,Spots,rot_start,rot_step, ...
        S,B,Ahkl,nrhkl,hkl_square,Energy,lambda,V,K1,I0E,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det, ...
        dety00,detz00,P0y,P0z,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        simap_data_flag,strcat(OutputFolder,'/grains_forfit'),rot_angles_simu);
if length(DS_new.SeedID)==length(DS_fit.SeedID)
    SpotsPair_all=SpotsPair;
end
SpotsPair_all0=SpotsPair_all;
SpotsPair0=SpotsPair;
toc

% make spots pair unique and select the good ones for fitting
SpotsPair=make_unique_spot_pair(SpotsPair0);
SpotsPair_all=make_unique_spot_pair(SpotsPair_all0);
%     SpotsPair=SpotsPair(SpotsPair(:,20)<25,:);

% col 13-14 forward simulated position; col 16-17 experimental position [pixel]
[ErrMean0,Err0,dis_y0,dis_z0]=dis_calc_spotspair([Lsam2sou Lsam2det dety00 detz00 tilt_x tilt_y tilt_z], ...
    SpotsPair,S,B,P0y,P0z,pixelysize,pixelzsize,dety0,detz0);
[ErrMean0_all,Err0_all,dis_y0_all,dis_z0_all]=dis_calc_spotspair([Lsam2sou Lsam2det dety00 detz00 tilt_x tilt_y tilt_z], ...
    SpotsPair_all,S,B,P0y,P0z,pixelysize,pixelzsize,dety0,detz0);
fprintf('%d spots are selected for fitting the geometry. \n',length(SpotsPair(:,1)));
fprintf('Average distance for selected grains is %.2f pixels. \n',ErrMean0);
fprintf('Average distance for all grains is %.2f pixels. \n',ErrMean0_all);

% add on Nov 25, 2021, random testing on different initial values
iter_fit=1;
for j=1:iter_fit
    if j==1
        x0(j,:)=[Lsam2sou Lsam2det dety00 detz00 tilt_x tilt_y tilt_z];
    else
        x0(j,1)=Lsam2sou-1+rand(1)*2;
        x0(j,2)=Lsam2det-1+rand(1)*2;
        x0(j,3)=dety00-0.1+rand(1)*0.2;
        x0(j,4)=detz00-0.2+rand(1)*0.2;
        x0(j,5)=tilt_x-1+rand(1)*2;
        x0(j,6)=tilt_y-1+rand(1)*2;
        x0(j,7)=tilt_z-1+rand(1)*2;
    end
    [ParaFit(j,:),fval(j)]=L_shift_tilt_fitting(SpotsPair,S,B,P0y,P0z,pixelysize,pixelzsize, ...
        dety0,detz0,x0(j,1),x0(j,2),x0(j,3),x0(j,4),x0(j,5),x0(j,6),x0(j,7),'FitAllOnce');
    j
end
[value,ind]=min(fval);

ParaFit=ParaFit(ind,:);
fval=fval(ind);

% col 13-14 forward simulated position; col 16-17 experimental position [pixel]
[ErrMean_fit,Err_fit,dis_y_fit,dis_z_fit]=dis_calc_spotspair(ParaFit,SpotsPair,S,B,P0y,P0z,pixelysize,pixelzsize,dety0,detz0);
[ErrMean_all_fit,Err_all_fit,dis_y_all_fit,dis_z_all_fit]=dis_calc_spotspair(ParaFit, ...
    SpotsPair_all,S,B,P0y,P0z,pixelysize,pixelzsize,dety0,detz0);
fprintf('Average distances for selected grains before and after fitting are %.2f and %.2f pixels,respectively. \n', ...
    ErrMean0,ErrMean_fit);
fprintf('Average distances for all grains before and after fitting are %.2f and %.2f pixels,respectively. \n', ...
    ErrMean0_all,ErrMean_all_fit);

Forward_simu_spots_exp_after_geo_fit(DS_fit,ParaFit, ...
    Rsample,RecVolumePixel,tomo_scale,ExpTime,atomparam,proj,Spots,rot_start,rot_step, ...
    S,B,Ahkl,nrhkl,hkl_square,Energy,lambda,V,K1,I0E,thetamax,lambda_min,lambda_max, ...
    P0y,P0z,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
    simap_data_flag,OutputFolder,[rot_start:2*rot_step:rot_end-180]);

if true
    plot_dis;
%     slice0=slice_show(DS_new,round(DS_new.Dimension(3)/2),1,0,0)
end
