function [] = Dream3DWriter(DS,dfile)
%%
% The Dream3DWriter function writes a .dream3d file 
% Parameters
% dfile: string of path to the dream3d file
% DS: struct dataset containing the result from GrainMapping

% MTEX toolbox is required and its path should be added
% Modified by Haixing Fang on March 22, 2019

datatype = containers.Map({'float','uint32','string','int32'}, {['DataArray<float>',0],['DataArray<uint32_t>',0], ['StringDataArray',0],['DataArray<int32_t>',0]});
axdim = {@(x) sprintf(['x=%d',0],x); @(x,y) sprintf(['x=%d,y=%d',0],x,y);@(x,y,z) sprintf(['x=%d,y=%d,z=%d',0],x,y,z)};

dataSetName = 'GM3D';
%%calculate the euler angle from the Rodrigues vector.

euler = Euler(rodrigues2quat(vector3d(double(DS.RodVec'))));
% % if there is no mtex toolbox
% [yaw, pitch, roll]=rod2angle(DS.RodVec,'ZXZ');
% yaw=yaw*180/pi;
% if yaw<0
%     yaw=yaw+360;
% end
% pitch=pitch*180/pi;
% if pitch<0
%     pitch=pitch+360;
% end
% roll=roll*180/pi;
% if roll<0
%     roll=roll+360;
% end
% euler = [yaw pitch roll];

EulerAngles = zeros([3 size(DS.GIDvol)]);
for i = 1:length(DS.SeedID)
    lind = find(DS.GIDvol == i);
    [dx,dy,dz] = ind2sub(size(DS.GIDvol),lind);
    for j = 1:length(lind)
        EulerAngles(:,dx(j),dy(j),dz(j)) = euler(i,:);
    end
end

featureNUM = length(DS.SeedID);
%% Prepare the file
% Create the file and empty groups in the file

if exist(dfile, 'file')~=2
    fid = H5F.create(dfile);% file path for creating the dream3D file
else
    dfile=fullfile(pwd,'New.dream3d'); 
    fid = H5F.create(dfile);
end
plist = 'H5P_DEFAULT';
gid = H5G.create(fid,'DataContainerBundles',plist,plist,plist);
H5G.close(gid);
gid = H5G.create(fid,'Pipeline',plist,plist,plist);
H5G.close(gid);
H5F.close(fid);
h5writeatt(dfile,'/Pipeline','Number_Filters',int32(0));

% add file info to root /
h5writeatt(dfile,'/','DREAM3D Version',['6.2.327.4fc644',0]);
h5writeatt(dfile,'/','FileVersion',['7.0',0]);

%% Write the cell data
dx = DS.Dimension(1); 
dy = DS.Dimension(2); 
dz = DS.Dimension(3);

% create and write completeness
comp_dream3d = permute(DS.CompVol, [4,1,2,3]);

h5create(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),size(comp_dream3d),'DataType','single');
h5write(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),single(comp_dream3d));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),'ObjectType',datatype('float'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Completeness',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

if isfield(DS,'Dismedian')
    % create and write Dismedian
    DisMed_dream3d = permute(DS.Dismedian, [4,1,2,3]);
    h5create(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),size(DisMed_dream3d),'DataType','single');
    h5write(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),single(DisMed_dream3d));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),'ComponentDimensions',uint64(1));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),'DataArrayVersion',int32(2));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),'ObjectType',datatype('float'));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Dismedian',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));
end
% write CellData attributes
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData',dataSetName),'AttributeMatrixType',uint32(3));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));



% create and write mask
mask_dream3d = permute(DS.Mask, [4,1,2,3]);

h5create(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),size(mask_dream3d),'DataType','int32');
h5write(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),int32(mask_dream3d));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),'ObjectType',datatype('int32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Mask',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

h5writeatt(dfile,sprintf('/DataContainers/%s/CellData',dataSetName),'AttributeMatrixType',uint32(3));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

% create and write phase
h5create(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),size(mask_dream3d),'DataType','int32');
h5write(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),ones(size(mask_dream3d),'int32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),'ObjectType',datatype('int32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/Phases',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

%create and write grainID
grainID_dream3d = permute(DS.GIDvol, [4,1,2,3]);

h5create(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),size(grainID_dream3d),'DataType','int32');
h5write(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),int32(grainID_dream3d));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),'ObjectType',datatype('int32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/GrainID',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

%create and write Euler angles data

h5create(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),size(EulerAngles),'DataType','single');
h5write(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),single(EulerAngles));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),'ComponentDimensions',uint64(3));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),'ObjectType',datatype('float'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/EulerAngles',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));

%create and write IPF data:
if isfield(DS,'IPF001')
    IPF001 = DS.IPF001;
    h5create(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),size(IPF001),'DataType','single');
    h5write(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),single(IPF001));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),'ComponentDimensions',uint64(3));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),'DataArrayVersion',int32(2));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),'ObjectType',datatype('float'));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),'Tuple Axis Dimensions',axdim{3}(dx,dy,dz));
    h5writeatt(dfile,sprintf('/DataContainers/%s/CellData/IPF001',dataSetName),'TupleDimensions',uint64([dx,dy,dz]));
end   
    

%% Create cell ensemble and write attributes
h5create(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),[1, 2],'DataType','uint32');
h5write(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),uint32([999, 1]));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),'ObjectType',datatype('uint32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),'Tuple Axis Dimensions',axdim{1}(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/CrystalStructures',dataSetName),'TupleDimensions',uint64(2));

h5create(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),[6,2],'DataType','single');
LatticeConstants = single([0,0,0,0,0,0; 2.8665, 2.8665, 2.8665, 90, 90, 90]);
h5write(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),LatticeConstants');
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),'ComponentDimensions',uint64(6));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),'ObjectType',datatype('float'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),'Tuple Axis Dimensions',axdim{1}(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/LatticeConstants',dataSetName),'TupleDimensions',uint64(2));

%% write Cellfeaturedata

h5create(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),size(DS.nVox),'DataType','uint32');
% h5create(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),size(DS.nVox),'DataType','uint32','ChunkSize',size(DS.nVox));% Dec 23, 2020
h5write(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),uint32(DS.nVox));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),'ObjectType',datatype('int32'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),'Tuple Axis Dimensions',axdim{1}(featureNUM));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellFeatureData/nVox',dataSetName),'TupleDimensions',uint64(featureNUM));

h5writeatt(dfile,'/DataContainers/GM3D/CellFeatureData','AttributeMatrixType',uint32(7));
h5writeatt(dfile,'/DataContainers/GM3D/CellFeatureData','TupleDimensions',uint64(featureNUM));


%% Write cellensemble attributes
MaterialName = {'Invalid phase';'Ferrite'};
plist = 'H5P_DEFAULT';
fid = H5F.open(dfile,'H5F_ACC_RDWR',plist);
% Set variable length string type
VLstr_type = H5T.copy('H5T_C_S1');
H5T.set_size(VLstr_type,'H5T_VARIABLE');
dspace = H5S.create_simple(1,numel(MaterialName),numel(MaterialName));

% Create dataset
dset = H5D.create(fid,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),VLstr_type,dspace,plist);
% Write data
H5D.write(dset,VLstr_type,'H5S_ALL','H5S_ALL','H5P_DEFAULT',MaterialName);
% Close file & resources
H5P.close(plist);
H5T.close(VLstr_type);
H5S.close(dspace);
H5D.close(dset);
H5F.close(fid);
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),'ComponentDimensions',uint64(1));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),'DataArrayVersion',int32(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),'ObjectType',datatype('string'));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),'Tuple Axis Dimensions',axdim{1}(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData/MaterialName',dataSetName),'TupleDimensions',uint64(2));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData',dataSetName),'AttributeMatrixType',uint32(11));
h5writeatt(dfile,sprintf('/DataContainers/%s/CellEnsembleData',dataSetName),'TupleDimensions',uint64(2));

% create and write _simpl_geomtry
h5create(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/DIMENSIONS',dataSetName),3,'DataType','int64');
h5write(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/DIMENSIONS',dataSetName),int64([dx,dy,dz]));
h5create(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/ORIGIN',dataSetName),3,'DataType','single');
h5write(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/ORIGIN',dataSetName),single([0,0,0]));
h5create(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/SPACING',dataSetName),3,'DataType','single');
h5write(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY/SPACING',dataSetName),single(DS.VoxSize));
h5writeatt(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY',dataSetName),'GeometryName',['ImageGeometry',0]);
h5writeatt(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY',dataSetName),'GeometryType',uint32(11));
h5writeatt(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY',dataSetName),'GeometryTypeName',['ImageGeometry',0]);
h5writeatt(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY',dataSetName),'SpatialDimensionality',uint32(3));
h5writeatt(dfile,sprintf('/DataContainers/%s/_SIMPL_GEOMETRY',dataSetName),'UnitDimensionality',uint32(3));
