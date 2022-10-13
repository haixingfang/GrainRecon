% get parameters
% May 31, 2022

% geometry information
Lsam2sou=geo.Lsam2sou;
Lsam2det=geo.Lsam2det;
P0y=geo.P0y;
P0z=geo.P0z;
if isfield(geo,'RotAxisOffset')
    RotAxisOffset=geo.RotAxisOffset;
end
dety00=geo.dety00;
detz00=geo.detz00;
tilt_x=geo.tilt_x;
tilt_y=geo.tilt_y;
tilt_z=geo.tilt_z;
RotDet=geo.RotDet;

% detector information
pixelysize=det.pixelysize;
pixelzsize=det.pixelzsize;
dety0=det.dety0;
detz0=det.detz0;
detysize=det.detysize;
detzsize=det.detzsize;
BeamStopY=det.BeamStopY;
BeamStopZ=det.BeamStopZ;

% reconstruction parameters 
drop_off=rec.drop_off;
minComp=rec.minComp;
minEucDis=rec.minEucDis;
maxD=rec.maxD;
RecVolumePixel=rec.RecVolumePixel;
VoxSize=rec.VoxSize;  % [mm]
TrustComp=rec.TrustComp;
Ahkl=rec.Ahkl;
OR=rec.OR;
if isfield(rec,'maxDmedian')
    maxDmedian=rec.maxDmedian;
end

% sample information
FileFolder=sample.FileFolder;
OutputFolder=sample.OutputFolder;
simap_data_flag=sample.simap_data_flag;
sample_no=sample.sample_no;
fitted_geo_already=sample.fitted_geo_already;
rot_angles=sample.rot_angles;
tomoFile=sample.tomoFile;
if isfield(sample,'SpotsFile')
    SpotsFile=sample.SpotsFile;
end
if ~exist(OutputFolder,'dir')
    mkdir(OutputFolder);
end

