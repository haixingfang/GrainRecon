% forward simulation
% sample: virtual_Fe_100um_6grains, AlCu_8wt_middle_thinned, SR-DCT measured on Sep 15, 2021
% Oct 6, 2021

clear all;
close all;

% % ground truth
FileFolder='./Examples/virtual_Fe_100um_6grains';
h5FileName='Grain100um_6grains_150_150_200_input.h5';
DS=readLabDCT(fullfile(FileFolder,h5FileName));

OutputFolder=FileFolder;
% load cell parameters;
input_fe;
B=FormB(cell);
V = cellvolume(cell); % [Angs^3]
simap_data_flag=1; % using simap data = 1; otherwise = 0
sample_no=4; % sample number
SR_flag=0; % data from synchrotron DCT ?
fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no

% load experimental parameters: source, detector and distances
setup_exp;
L=Lsam2sou+Lsam2det;
RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
RotDet=RotX*RotY*RotZ;

% generate hkl for indexing
hklnumber = 4; % choose how many hkl families to use for simulation
setup_hkl;

% other parameters
if simap_data_flag==1
    S=[1 0 0;0 1 0;0 0 1];
else
    S=[1 0 0;0 -1 0;0 0 1];
end
SampleVolumeDim=[0.15 0.15 0.2]; % [mm]
Rsample=sqrt(SampleVolumeDim(1)^2+SampleVolumeDim(2)^2)/2; % [mm]

RecVolumePixel(:,1)=[1 1 1]';
RecVolumePixel(:,2)=size(DS.GIDvol)';
DIM=size(DS.GIDvol);

% % grain verification and pairing spots
grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]

% to create simulated images in simap geometry
rot_start=0;
rot_step=3;
rot_end=360;
rot_angles=rot_start:rot_step:rot_end;
% rot_angles_simu=[117];
rot_angles_simu=rot_angles;
simu_grainno=[find(grainsize>0)]';
% simu_grainno=[1 2 3 4];
fprintf('Start computing simulated images ...\n')
tic
[SpotNr_simu,GrainIndex_all,SubGrain,rot_angles_simulated,SmallGrID]=Forward_simu(DS,Rsample,RecVolumePixel, ...
    DIM,ExpTime,atomparam,rot_start,rot_step,rot_end,S,B,Ahkl,nrhkl,hkl_square, ...
    Energy,lambda,V,K1,I0E,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
    simap_data_flag,OutputFolder,SR_flag,rot_angles_simu,simu_grainno);
t = toc;
fprintf('Done! Total computation time is %.1f s.\n', t);
