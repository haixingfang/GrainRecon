% get tomo slices and create a h5 file for the tomo volume
% % input: tomo slices from Xact;
%          xml file from Xact tomo reconstruction
% August 17, 2021

clear all;
close all;

%%%%% acquisition file
% xmlFolder='D:\ExpData\2021_07_22_AlCu8wt\tomo-0723-no-pinhole';
% xmlFolder='D:\ExpData\2021_09_30_DCT\AlCu_8wt_middle_thinned_Lab\real_scan\abs_tomo';
% xmlFolder='D:\ExpData\2022_01_19_DCT\AlCu_8wt_middle_thinned_micro_source_revise\tomo_no_pinhole';
% xmlFolder='D:\ExpData\2022_01_19_DCT\Al8Cu_pin_580C_1h\tomo_no_pinhole';
xmlFolder='D:\ExpData\2022_01_19_DCT\Al8Cu_pin_580C_1h\tomo_no_pinhole';
xmlName='unireconstruction.xml';
xmlfile=xml2struct(fullfile(xmlFolder,xmlName));
voxelsize=mean([str2num(xmlfile.unireconstruction.conebeam.volume.voxelSize.Attributes.X) ...
    str2num(xmlfile.unireconstruction.conebeam.volume.voxelSize.Attributes.Y) ...
    str2num(xmlfile.unireconstruction.conebeam.volume.voxelSize.Attributes.Z)])*1000; % [um]

%%%% get slices
% ImageFolder='D:\ExpData\2021_07_22_AlCu8wt\tomo-0723-no-pinhole\SlicesY';
% ImageFolder='D:\ExpData\2021_09_30_DCT\AlCu_8wt_middle_thinned_Lab\real_scan\abs_tomo\SlicesY';
% ImageFolder='D:\ExpData\2022_01_19_DCT\AlCu_8wt_middle_thinned_micro_source_revise\tomo_no_pinhole\SlicesY';
% ImageFolder='D:\ExpData\2022_01_19_DCT\Al8Cu_pin_580C_1h\tomo_no_pinhole\SlicesY';
ImageFolder='D:\ExpData\2022_01_19_DCT\Al8Cu_pin_580C_1h\tomo_no_pinhole\SlicesY';
name1='slice';
name3='.tif';
files = dir([ImageFolder '\*',name3]);
MinGrainSize=100;

for i=1:numel(files)
    imdata(:,:,i)=imread(fullfile(ImageFolder,files(i).name));
    [im_size(1),im_size(2),~]=size(imdata);
end
figure;view0=orthosliceViewer(imdata);
ThresholdValue=22000; % by default 25500
im_bin=imdata>ThresholdValue;
im_bin=bwareaopen(im_bin,MinGrainSize);
im_label=bwlabeln(im_bin,18);
im_msr=regionprops3(im_label,'all');
figure;view1=orthosliceViewer(im_label,'Colormap',parula(256),'DisplayRange',[0 10]);

objectID=find(im_msr.Volume>(500/voxelsize)^3*0.3); % sample volume
% to select only the cubes
im_bin=zeros(size(im_bin));
for j=1:length(objectID)
    ID=objectID(j);
    im_select=im_label==ID;
    im_bin=im_bin+im_select;
    j
end
im_bin=im_bin>0;
im_label=bwlabeln(im_bin,18); % label again
im_msr=regionprops3(im_label,'all'); % measure again
figure;view2=orthosliceViewer(im_label,'Colormap',parula(256),'DisplayRange',[0 10]);

%%%% write tomo h5 file
% dimensions of x, y, z, usual coordinate
% OutputFolder='D:\ExpData_process\CCD_2021_07_22_AlCu8wt';
% OutputFolder='D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned';
% OutputFolder='./AlCu_8wt_middle_thinned_micro_source';
OutputFolder='./AlCu_8wt_580C_1h_micro_source';
if ~exist(OutputFolder, 'dir')
   mkdir(OutputFolder);
end

grains_out=im_label;
dim=size(grains_out);

% from the reconstruction file
Center0=[str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.X) ...
    str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Y) ...
    str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Z)]; % [mm]
Center_crop=[str2num(xmlfile.unireconstruction.conebeam.volume.offset.Attributes.X) ...
    str2num(xmlfile.unireconstruction.conebeam.volume.offset.Attributes.Y) ...
    str2num(xmlfile.unireconstruction.conebeam.volume.offset.Attributes.Z)]; % [mm]
Center(1)=Center_crop(3)-Center0(3);
Center(2)=Center_crop(1)-Center0(1);
Center(3)=Center_crop(2)-Center0(2);

% write tomo h5
% FileName='tomo_AlCu8wt_no_pinhole_0723.h5';
% FileName='tomo_2021_09_30_AlCu_8wt_middle_thinned.h5';
% FileName='tomo_2022_01_19_AlCu_8wt_middle_thinned_microS.h5';
FileName='tomo_2022_01_20_AlCu_8wt_580C_1h_microS.h5';
ProjectName=FileName(1:end-3);
tomo=write_h5tomo_from_im_label(grains_out,Center,[],[],'Al fcc', ...
    voxelsize,OutputFolder,FileName,ProjectName);
Prefix_Name=FileName(1:end-3);
save(fullfile(OutputFolder,[Prefix_Name '.mat']));



% % for simulated vol, consider use the follows
% im=[];
% for i=1:size(tomo.Mask,3)
%     im(:,:,i)=imrotate(tomo.Mask(:,:,i),90);
% end
% grains_out=im;
% FileName='tomo_virtual_Fe_100um_6grains_new.h5';
% ProjectName=FileName(1:end-3);
% tomo_new=write_h5tomo_from_im_label(grains_out,Center,[],[],'Al fcc', ...
%     voxelsize,OutputFolder,FileName,ProjectName);

%{
% for fast checking the sample position
DS.Center=Center;
DS.Dimension=dim;
DS.VoxSize=[voxelsize voxelsize voxelsize]/1000;
GrainID=uint16(grains_out);
for i=1:length(im_msr.Volume)
    [x,y,z] = ind2sub(size(GrainID),find(GrainID == i));
    X = mean(x);
    Y = mean(y);
    Z = mean(z);
    DS.Coord(i,:) = [X,Y,Z];
end
center_shift(1)= DS.Center(1); % x: along the beam [mm]   (lab-z)
center_shift(2)= DS.Center(2); % y: perpendicular to the beam [mm] (lab-x)
center_shift(3)= DS.Center(3); % z: sample height direction, irrelevant [mm] (lab-y)
sampos(:,1) = (DS.Coord(:,1)-DS.Dimension(1)/2).*DS.VoxSize(1) ...
    +center_shift(1); % coordinates of grain centroids [mm]
sampos(:,2) = (DS.Coord(:,2)-DS.Dimension(2)/2).*DS.VoxSize(2) ...
    +center_shift(2); % coordinates of grain centroids [mm]
sampos(:,3) = (DS.Coord(:,3)-DS.Dimension(3)/2).*DS.VoxSize(3)+center_shift(3);    
sampos(:,1)=-sampos(:,1);  % corrected on July 6th, 2021
sampos(:,2)=-sampos(:,2);  % corrected on July 6th, 2021
%}
