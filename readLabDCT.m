function resultFile = readLabDCT(filepath)

%%
% The readLabDCT function reads the result file from GrainMapper3D and
% write into a struct file.
% Parameters
% filepath: string of path to the result file

% modified by Haixing Fang on March 22, 2019

fid = H5F.open(filepath);
fapl = H5F.get_access_plist(fid);

DataInfo = h5info(filepath);
for i=1:length(DataInfo.Groups(2).Datasets) 
    PropertyName{i} = DataInfo.Groups(2).Datasets(i).Name;
end
str0=strjoin(PropertyName,',')
for i=1:length(DataInfo.Groups(2).Groups.Datasets) 
    PropertyName{i} = DataInfo.Groups(2).Groups.Datasets(i).Name;
end
str=strjoin(PropertyName,',')

resultFile.Center = h5read(filepath,'/LabDCT/Center');
resultFile.VoxSize = h5read(filepath,'/LabDCT/Spacing');
resultFile.GIDvol = h5read(filepath,'/LabDCT/Data/GrainId');
resultFile.CompVol = h5read(filepath,'/LabDCT/Data/Completeness');
% corrected on Feb 4, 2022
if ~isempty(strfind(str, 'Mask'))
    resultFile.Mask = h5read(filepath,'/LabDCT/Data/Mask'); 
else
    resultFile.Mask = h5read(filepath,'/LabDCT/Data/PhaseId');
end
resultFile.Dimension = size(resultFile.GIDvol);

if ~isempty(strfind(str0, 'RecVolumePixel'))
    resultFile.RecVolumePixel = h5read(filepath,'/LabDCT/RecVolumePixel');
else
    fprintf('/LabDCT/RecVolumePixel property does not exist\n');
end
if ~isempty(strfind(str0, 'tomo_scale_dim'))
    resultFile.tomo_scale_dim = h5read(filepath,'/LabDCT/tomo_scale_dim');
else
    fprintf('/LabDCT/tomo_scale_dim property does not exist\n');
end

if ~isempty(strfind(str, 'Rodrigues'))
    resultFile.RodVec3D = h5read(filepath,'/LabDCT/Data/Rodrigues');
else
    fprintf('/LabDCT/Data/Rodrigues property does not exist\n');
end
if ~isempty(strfind(str, 'EulerZXZ'))
    resultFile.EulerAngle = h5read(filepath,'/LabDCT/Data/EulerZXZ');
else
    fprintf('/LabDCT/Data/EulerZXZ property does not exist\n');
end
if ~isempty(strfind(str, 'IPF001'))
    resultFile.IPF001 = h5read(filepath,'/LabDCT/Data/IPF001');
else
    fprintf('/LabDCT/Data/IPF001 property does not exist\n');
end
if ~isempty(strfind(str, 'Dismedian'))
    resultFile.Dismedian = h5read(filepath,'/LabDCT/Data/Dismedian');
else
    fprintf('/LabDCT/Data/Dismedian property does not exist\n');
end

resultFile.SeedID = 1:max(max(max(resultFile.GIDvol)));


for i = resultFile.SeedID
%     resultFile.SeedComp(i,1) = mean(resultFile.CompVol(resultFile.GIDvol == i)); % completeness for each grain
    resultFile.SeedComp(i,1) = nanmean(resultFile.CompVol(resultFile.GIDvol == i)); % completeness for each grain
end


if ~isempty(strfind(str, 'Rodrigues'))
    for i = resultFile.SeedID
        ex = resultFile.RodVec3D(:,resultFile.GIDvol == i);
        if isempty(ex)
            resultFile.RodVec(i,:) = [0 0 0];
        else
            resultFile.RodVec(i,:) = ex(:,1)'; % each row represents Rodrigues vector for each grain
        end
    end
end
if ~isempty(strfind(str, 'EulerZXZ'))
    for i = resultFile.SeedID
        ex = resultFile.EulerAngle(:,resultFile.GIDvol == i);
        if isempty(ex)
            resultFile.EulerZXZ(i,:) = [0 0 0];
        else
            resultFile.EulerZXZ(i,:) = ex(:,1)'; % each row represents ZXZ Euler angles for each grain
            if sum(ex(:,1)'>2*pi)==0
                resultFile.EulerZXZ(i,:)=resultFile.EulerZXZ(i,:)*180/pi; % [degrees]
            end
        end
    end
end


for i = resultFile.SeedID
    resultFile.nVox(i,1) = length(find(resultFile.GIDvol == i)); % number of voxels for each grain
end

for i = resultFile.SeedID
    [x,y,z] = ind2sub(size(resultFile.GIDvol),find(resultFile.GIDvol == i));
    X = mean(x);
    Y = mean(y);
    Z = mean(z);
    resultFile.Coord(i,:) = [X,Y,Z]; % each row represents coodinate for each grain
end



