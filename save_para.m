% save parameters
% Nov 19, 2021
function save_para(Lsam2sou,Lsam2det,P0y,P0z,dety00,detz00,tilt_x,tilt_y,tilt_z,RotDet, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
    drop_off,maxD,TrustComp,minComp,minEucDis,Ahkl,OR,VoxSize,RecVolumePixel, ...
    FileFolder,OutputFolder,simap_data_flag,sample_no,fitted_geo_already,rot_angles, ...
    tomo_scale,h5Folder_tomo,h5FileName_tomo,SpotsFile,maxDmedian)

% geometry information
geo.Lsam2sou=Lsam2sou;
geo.Lsam2det=Lsam2det;
geo.P0y=P0y;
geo.P0z=P0z;
geo.dety00=dety00;
geo.detz00=detz00;
geo.tilt_x=tilt_x;
geo.tilt_y=tilt_y;
geo.tilt_z=tilt_z;
geo.RotDet=RotDet;

% detector information
det.pixelysize=pixelysize;
det.pixelzsize=pixelzsize;
det.dety0=dety0;
det.detz0=detz0;
det.detysize=detysize;
det.detzsize=detzsize;
det.BeamStopY=BeamStopY;
det.BeamStopZ=BeamStopZ;

% reconstruction parameters 
rec.drop_off=drop_off;
rec.minComp=minComp;
rec.minEucDis=minEucDis;
rec.maxD=maxD;
rec.RecVolumePixel=RecVolumePixel;
rec.VoxSize=VoxSize;  % [mm]
rec.TrustComp=TrustComp;
rec.Ahkl=Ahkl;
rec.OR=OR;
if nargin>37
    rec.maxDmedian=maxDmedian;
end

% sample information
sample.FileFolder=FileFolder;
sample.OutputFolder=OutputFolder;
sample.simap_data_flag=simap_data_flag;
sample.sample_no=sample_no;
sample.fitted_geo_already=fitted_geo_already;
sample.rot_angles=rot_angles;
sample.center=tomo_scale.Center; % [mm]
sample.VoxSize=tomo_scale.VoxSize; % [mm]
sample.tomoFile=fullfile(h5Folder_tomo,h5FileName_tomo);
if nargin>36
    sample.SpotsFile=fullfile(FileFolder,SpotsFile);
end

if ~exist(OutputFolder,'dir')
    mkdir(OutputFolder);
end
save(fullfile(OutputFolder,'para.mat'),'geo','det','rec','sample');




