% transform im_label to h5 file for absorption tomographic volume
% export dream3d and xdmf file for 3d view (optional)
% July 8, 2021
function tomo=write_h5tomo_from_im_label(grains_out,Center,CenterShift,PhaseNo,PhaseName, ...
    voxelsize,OutputFolder,FileName,ProjectName)

% dimensions of x, y, z, usual coordinate
% grains_out=im_label;
dim=size(grains_out); 
if isempty(CenterShift)
    CenterShift=[0 0 0];
end
if isempty(PhaseNo)
    PhaseNo=1;
end

pixelsize=voxelsize/1000; % [mm]
Spacing=[pixelsize pixelsize pixelsize]';
VirtualShift=[0 0 0]';
Extent=dim'*pixelsize;

% initialize
Mask=zeros(dim);
PhaseId=zeros(dim);
grains_out=grains_out>0;
clear indices;
[indices(:,1),indices(:,2),indices(:,3)]=ind2sub(size(grains_out),find(grains_out==1));
if ~isempty(indices)
    for j=1:length(indices(:,1))
        PhaseId(indices(j,1),indices(j,2),indices(j,3))=PhaseNo;
        Mask(indices(j,1),indices(j,2),indices(j,3))=1;
    end
end
hdf5write(fullfile(OutputFolder,FileName), '/AbsorptionCT/Date', datestr(now,'dd/mm/yy'), ...
    '/AbsorptionCT/Center', double(Center), ...
    '/AbsorptionCT/CenterShift', double(CenterShift), ...
    '/AbsorptionCT/Spacing', double(Spacing), ...
    '/AbsorptionCT/Data/PhaseId', uint8(PhaseId), ...
    '/AbsorptionCT/Data/Mask', double(Mask), ...
    '/AbsorptionCT/VirtualShift', double(VirtualShift), ...
    '/AbsorptionCT/Extent', double(Extent), ...
    '/PhaseInfo/Phase01/Name', PhaseName, ...
    '/ProjectInfo/ProjectFile', ProjectName);
% H5F.close(FileName);
tomo.Center = Center;
tomo.CenterShift = CenterShift;
tomo.VirtualShift = VirtualShift;
tomo.Dimension = dim;
tomo.VoxSize = Spacing;
tomo.Mask = Mask;
tomo.PhaseId = PhaseId;
tomo.PhaseName = PhaseName;


% FileFolder=OutputFolder;
% Prefix_Name=FileName(1:end-3);
% Dream3D_FileName=[Prefix_Name '.dream3d']; % create a new name for dream3D
% dfile=fullfile(FileFolder,Dream3D_FileName); % file path for creating the dream3D file
% Dream3DWriter_tomo(tomo,dfile);






