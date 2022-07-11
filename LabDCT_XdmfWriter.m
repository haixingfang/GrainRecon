% write a Xdmf file for visualization by Paraview from a dream3D file
% written by Haixing Fang March 22, 2019

function LabDCT_XdmfWrite(FileFolder,Dream3D_FileName,Xdmf_FileName)

% % example for input:
% % FileFolder='D:\Users\hfang\AL_APSESRF_analysis\Dss-14_Dds14_Exp600s_Feb24 2017'; % Avizo PC
% FileFolder='M:\MetaData\Al_Dss-14_Dsd14_exp600s_Feb24_2017'; % HP laptop
% Dream3D_FileName='3Dvolume_4grains.dream3d';
% Xdmf_FileName='3Dvolume_4grains.xdmf';

Dream3D_fullFileName=fullfile(FileFolder,Dream3D_FileName);
Xdmf_fullFileName=fullfile(FileFolder,Xdmf_FileName);

GIDvol = h5read(Dream3D_fullFileName,'/DataContainers/GM3D/CellData/GrainID');
dimensions = flip(size(GIDvol)); % dimensions of the dataset
MeshDimensions = dimensions(1:3)+1; % dimensions of the mesh
ScalarDimensions = dimensions; % dimensions of scalar data: completeness, GrainID, Mask, Phases
VectorDimensions(1:3) = dimensions(1:3); % dimensions of vector data
VectorDimensions(4) = 3; % dimensions of vector data: Euler Angles, IPF001

voxelsize = h5read(Dream3D_fullFileName,'/DataContainers/GM3D/_SIMPL_GEOMETRY/SPACING'); % [mm]

MeshDimensionsStr=sprintf('Dimensions="%d %d %d "',dimensions(1),dimensions(2),dimensions(3));
ScalarDimensionsStr=sprintf('Dimensions="%d %d %d %d"',ScalarDimensions);
VectorDimensionsStr=sprintf('Dimensions="%d %d %d %d"',VectorDimensions);

VoxelStr=sprintf('%.3f %.3f %.3f',voxelsize(1),voxelsize(2),voxelsize(3)); % [mm]

fileID=fopen(fullfile(FileFolder,Xdmf_FileName),'wt');
    fprintf(fileID,'<?xml version="1.0"?>\n');
    fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n');
    fprintf(fileID,'<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n');
    fprintf(fileID,' <Domain>\n');
    fprintf(fileID,'  <!-- *************** START OF GM3D *************** -->\n');
    fprintf(fileID,'  <Grid Name="GM3D" GridType="Uniform">');
    fprintf(fileID,['   <Topology TopologyType="3DCoRectMesh" ', MeshDimensionsStr, '></Topology>\n']);
    fprintf(fileID,'    <Geometry Type="ORIGIN_DXDYDZ">\n');
    fprintf(fileID,'     <!-- Origin  Z, Y, X -->\n');
    fprintf(fileID,'     <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>');
    fprintf(fileID,'     <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->');
    fprintf(fileID,['     <DataItem Format="XML" Dimensions="3">',VoxelStr,'</DataItem>']);
    fprintf(fileID,'    </Geometry>\n');
  
    fprintf(fileID,'    <Attribute Name="Completeness" AttributeType="Scalar" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',ScalarDimensionsStr,' NumberType="Float" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/Completeness','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'    <Attribute Name="EulerAngles" AttributeType="Vector" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',VectorDimensionsStr,' NumberType="Float" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/EulerAngles','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'    <Attribute Name="GrainID" AttributeType="Scalar" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',ScalarDimensionsStr,' NumberType="Int" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/GrainID','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'    <Attribute Name="IPF001" AttributeType="Vector" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',VectorDimensionsStr,' NumberType="Float" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/IPF001','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'    <Attribute Name="Mask" AttributeType="Scalar" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',ScalarDimensionsStr,' NumberType="Int" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/Mask','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'    <Attribute Name="Phases" AttributeType="Scalar" Center="Cell">\n');
    fprintf(fileID,['      <DataItem Format="HDF" ',ScalarDimensionsStr,' NumberType="Int" Precision="4" >', ...
        Dream3D_FileName,':','/DataContainers/GM3D/CellData/Phases','</DataItem>\n']);    
    fprintf(fileID,'    </Attribute>\n');
    
    fprintf(fileID,'  </Grid>\n');
    fprintf(fileID,'  <!-- *************** END OF GM3D *************** -->\n');
    fprintf(fileID,' </Domain>\n');
    fprintf(fileID,'</Xdmf>\n');
fclose(fileID);



