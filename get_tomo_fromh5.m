% Haixing Fang on July 8, 2021
function resultFile = get_tomo_fromh5(filepath,dataflag)
% dataflag=1: simap
% dataflag=0: Zeiss Xradia format

if dataflag==1
    fid = H5F.open(filepath);
    fapl = H5F.get_access_plist(fid);
    DataInfo = h5info(filepath);

    resultFile.Center = h5read(filepath,'/AbsorptionCT/Center');
    resultFile.CenterShift = h5read(filepath,'/AbsorptionCT/CenterShift');
    resultFile.VirtualShift = h5read(filepath,'/AbsorptionCT/VirtualShift');
    resultFile.VoxSize = h5read(filepath,'/AbsorptionCT/Spacing');
    resultFile.PhaseId = h5read(filepath,'/AbsorptionCT/Data/PhaseId');
    resultFile.Mask = h5read(filepath,'/AbsorptionCT/Data/Mask');
    resultFile.Dimension = size(resultFile.Mask);
    resultFile.PhaseName = h5read(filepath,'/PhaseInfo/Phase01/Name');
else
    %%% for DTU LabDCT format
    fid = H5F.open(filepath);
    fapl = H5F.get_access_plist(fid);
    DataInfo = h5info(filepath);

    resultFile.Center = h5read(filepath,'/LabDCT/Center');
    resultFile.CenterShift = [0 0 0];
    resultFile.VirtualShift = [0 0 0];
    resultFile.VoxSize = h5read(filepath,'/LabDCT/Spacing');
    resultFile.PhaseId = h5read(filepath,'/LabDCT/Data/PhaseId');
    resultFile.Mask = h5read(filepath,'/LabDCT/Data/PhaseId');
    resultFile.Dimension = size(resultFile.Mask);
    resultFile.PhaseName = h5read(filepath,'/PhaseInfo/Phase01/Name');
end
