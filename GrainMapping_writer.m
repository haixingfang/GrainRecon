% GrainMapping_writer
% write h5, dream3d and xdmf files from DS
% July 30, 2021
function [DS_new]=GrainMapping_writer(DS,Center,CenterShift,PhaseNo,PhaseName, ...
    voxelsize,RecVolumePixel,tomo_scale_dim,FileFolder,fname_prefix,ProjectName,regenerate_IPF)

if ~exist(FileFolder, 'dir')
   mkdir(FileFolder);
end

% dimensions of x, y, z, usual coordinate
dim=size(DS.Completeness); 
if isempty(CenterShift)
    CenterShift=[0 0 0];
end
if isempty(PhaseNo)
    PhaseNo=1;
end

if voxelsize<0.1
    pixelsize=voxelsize; % [mm]
else
    pixelsize=voxelsize/1000; % [mm]
end
Spacing=[pixelsize pixelsize pixelsize]';
VirtualShift=[0 0 0]';
Extend=dim'*pixelsize;

% initialize
PhaseId=DS.PhaseId*PhaseNo;
Mask=DS.Mask;
Completeness=DS.Completeness;
GrainId=DS.GrainId;
Rodrigues=DS.Rodrigues;
EulerZXZ=DS.EulerZXZ;
VisitFlag=DS.VisitFlag;
IPF001=DS.IPF001;
if isfield(DS,'Dismedian')
    Dismedian=DS.Dismedian;
end

if regenerate_IPF==1
    load_mtex;
    %%% change the default settings
    setMTEXpref('xAxisDirection','east'); % default: 'north'
    setMTEXpref('zAxisDirection','outOfPlane'); % same as default
    setMTEXpref('bAxisDirection','north'); % default: 'east'
    setMTEXpref('aAxisDirection',''); % same as default
    setMTEXpref('FontSize',44); % default: 15
    % define symmetries
    cs = crystalSymmetry('cubic');
    % ss = specimenSymmetry('orthorhombic');
    cK = ipfHSVKey(cs);
    cK.inversePoleFigureDirection=zvector;
    
    % faster
    SeedID = 1:max(max(max(DS.GrainId)));
    for i=1:length(SeedID)
        ind=find(DS.GrainId==i);
        if ~isempty(ind)
            ex=DS.EulerZXZ(:,ind);
            ex=ex(:,1)';       
            rot1=rotation('Euler',ex(1)*degree,ex(2)*degree,ex(3)*degree);
            o1=orientation(rot1,cs);
            GrainColor1 = cK.orientation2color(o1);
            IPF001(:,ind)=repmat(GrainColor1',1,length(ind));
        else
            sprintf('no seedID %d ...',i)
        end
%         i
        if mod(i,50)==0 || i==length(SeedID)
            sprintf('%d grains have been registered.', i)
        end
    end
    
%     % too slow
%     for i=1:length(DS.Completeness(:,1,1))
%         if all(nonzeros(DS.Completeness(i,:,:)))
%             for j=1:length(DS.Completeness(1,:,1))
%                 for k=1:length(DS.Completeness(1,1,:))
%                     if DS.Completeness(i,j,k)>0
%                         rot1=rotation('Euler',DS.EulerZXZ(1,i,j,k)*degree, ...
%                             DS.EulerZXZ(2,i,j,k)*degree,DS.EulerZXZ(3,i,j,k)*degree);
%     %                     o1=orientation(rot1,cs,ss);
%                         o1=orientation(rot1,cs);
%                         GrainColor1 = cK.orientation2color(o1);
%                         IPF001(:,i,j,k)=GrainColor1;
%                     else
%                         IPF001(:,i,j,k)=[NaN NaN NaN];
%                     end
%                 end
%             end
%         end
%         i
%     end
    
    DS.IPF001=IPF001;
end

h5FileName=[fname_prefix '.h5'];
if isfield(DS,'Dismedian')
    hdf5write(fullfile(FileFolder,h5FileName), '/AbsorptionCT/Date', datestr(now,'dd/mm/yy'), ...
        '/LabDCT/Center', double(Center), ...
        '/LabDCT/RecVolumePixel', double(RecVolumePixel'), ...
        '/LabDCT/tomo_scale_dim', tomo_scale_dim, ...
        '/LabDCT/Spacing', double(Spacing), ...
        '/LabDCT/Data/Mask', uint8(Mask), ...
        '/LabDCT/Data/GrainId', uint16(GrainId), ...
        '/LabDCT/Data/IPF001', double(IPF001), ...
        '/LabDCT/Data/Completeness', double(Completeness), ...
        '/LabDCT/Data/Dismedian', double(Dismedian), ...
        '/LabDCT/Data/PhaseId', uint8(PhaseId), ...
        '/LabDCT/Data/Rodrigues', double(Rodrigues), ...
        '/LabDCT/Data/EulerZXZ', double(EulerZXZ), ...
        '/LabDCT/Data/VisitFlag', double(VisitFlag), ...
        '/LabDCT/VirtualShift', double(VirtualShift), ...
        '/LabDCT/Extend', double(Extend), ...
        '/PhaseInfo/Phase01/Name', PhaseName, ...
        '/ProjectInfo/ProjectFile', ProjectName);
else
    hdf5write(fullfile(FileFolder,h5FileName), '/AbsorptionCT/Date', datestr(now,'dd/mm/yy'), ...
        '/LabDCT/Center', double(Center), ...
        '/LabDCT/RecVolumePixel', double(RecVolumePixel'), ...
        '/LabDCT/tomo_scale_dim', tomo_scale_dim, ...
        '/LabDCT/Spacing', double(Spacing), ...
        '/LabDCT/Data/Mask', uint8(Mask), ...
        '/LabDCT/Data/GrainId', uint16(GrainId), ...
        '/LabDCT/Data/IPF001', double(IPF001), ...
        '/LabDCT/Data/Completeness', double(Completeness), ...
        '/LabDCT/Data/PhaseId', uint8(PhaseId), ...
        '/LabDCT/Data/Rodrigues', double(Rodrigues), ...
        '/LabDCT/Data/EulerZXZ', double(EulerZXZ), ...
        '/LabDCT/Data/VisitFlag', double(VisitFlag), ...
        '/LabDCT/VirtualShift', double(VirtualShift), ...
        '/LabDCT/Extend', double(Extend), ...
        '/PhaseInfo/Phase01/Name', PhaseName, ...
        '/ProjectInfo/ProjectFile', ProjectName);
        % H5F.close(FileName);
end

DS = readLabDCT(fullfile(FileFolder,h5FileName));

Dream3D_FileName=[fname_prefix '.dream3d'];
dfile=fullfile(FileFolder,Dream3D_FileName); % file path for creating the dream3D file
if exist(dfile,'file')==0
    Dream3DWriter(DS,dfile);
else
    delete(dfile);
    Dream3DWriter(DS,dfile);
end
Xdmf_FileName=[Dream3D_FileName(1:end-8) '.xdmf']; % create a new name for xdmf file
if isfield(DS,'Dismedian')
    LabDCT_XdmfWriter_v2(FileFolder,Dream3D_FileName,Xdmf_FileName);
else
    LabDCT_XdmfWriter(FileFolder,Dream3D_FileName,Xdmf_FileName);
end
DS_new=DS;

