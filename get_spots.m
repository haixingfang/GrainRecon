% get Spot information from binarized images
% require: spot segmented images
% Aug 16, 2021
% function [proj,proj_bin,Spots,SpotsTable,rot_start,rot_end,rot_step,ImageNr]=SpotSegmentation(ImageFolder,name1,name3, ...
%     BeamStopY,BeamStopZ,Nrolling,BackGround,sigma,percent,MinSpotSize,write_im,parallelflag, ...
%     CombineFlag,CombineFlag_threshold)

clear all;
% ImageFolder='D:\ExpData_process\CCD_2021_07_22_AlCu8wt\DCT_1_6_all_proj_spot';
% OutputFileName='Spots_AlCu8wt_0722_121projs.mat';
% ImageFolder='D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\proj_spot';
% ImageFolder_label='D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\proj_label';
% ImageFolder='D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\proj_spot_Nov10';     % update on Nov 10
% ImageFolder_label='D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\proj_label_Nov10'; % update on Nov 10
% OutputFileName='Spots_AlCu_8wt_middle_thinned_0930.mat';

% ImageFolder='D:\ExpData_process\CCD_2022_01_19_AlCu_8wt_middle_thinned_micro_source\proj_spot';     % Jan 21, 2022
% ImageFolder_label='D:\ExpData_process\CCD_2022_01_19_AlCu_8wt_middle_thinned_micro_source\proj_label'; % Jan 21, 2022
% OutputFileName='Spots_AlCu_8wt_middle_thinned_micro_source.mat';

ImageFolder='D:\ExpData_process\CCD_2022_01_20_AlCu_8wt_pin_580C_1h\proj_spot';     % Jan 21, 2022
ImageFolder_label='D:\ExpData_process\CCD_2022_01_20_AlCu_8wt_pin_580C_1h\proj_label'; % Jan 21, 2022
OutputFileName='Spots_AlCu_8wt_580C_1h_micro_source.mat';

write_im_label=1;
name1='proj';
name3='.tif';
MinSpotSize=3; % 2-10 [pixels]

if exist(ImageFolder, 'dir')~=7
    Message = sprintf('Error: The following folder does not exist:\n%s', ImageFolder)
    %  uiwait(warndlg(Message));
    return;
end
if ~exist(ImageFolder_label,'dir')
    mkdir(ImageFolder_label);
end
ImageInfo=dir([ImageFolder '/*' name3]);
ImageNr = numel(ImageInfo);
rot_start=0;
rot_end=360;
rot_step=(rot_end-rot_start)/(ImageNr-1);
rot_number=1;
for rot=rot_start:rot_step:rot_end
    ProjectNo=(rot-rot_start)/rot_step+1;
    rotation_angle(rot_number)=rot; % [deg]
    name2=num2str(ProjectNo,'%.4d');
    ImageName=[name1 name2 name3];
    ImageName = fullfile(ImageFolder, ImageName);
    proj{rot_number}=readim(ImageName);
    if strcmp(name3,'.tiff')
        proj{rot_number}=flipud(proj{rot_number});
    end
    if max(size(proj{rot_number}))==3
        Im=dip_array(proj{rot_number},'double');
        Im=reshape(Im(:,:,1),size(Im,1),size(Im,2));
        proj{rot_number}=dip_image(Im,'uint16');
    end
    sprintf('Reading projection no. %d / %d',rot_number,ImageNr)
    rot_number=rot_number+1;
end

for i=1:length(proj)
    Spots{i}=[];
    im=proj{i};
    if i==1
       ImageSize=size(proj{i});
    end
    im_spot_bin=im>0; % when proj{i} is segmented already
    
%     BackGround=double(median(median(im)))+0.1; % modified on Dec 7, 2020
%     sigma=0.5;
% %     percent=0.075; % [0 0.25] minimum intensity fraction of a blob to be considered as a spot
%     percent=0.1; % March 8, 2021
%     MinSpotSize=2; % [pixels]
%     im_spot_bin=LoG_segmentation_parallel(im,BackGround,sigma,percent,MinSpotSize); % when proj{i} needs LoG segmentation  
    
    im_spot_bin=dip_image(im_spot_bin,'bin');
    proj_bin{i}=im_spot_bin;
    imdata=proj{i};
    im_label=label(im_spot_bin,2,MinSpotSize*2,size(im,1)*size(im,2));
    if write_im_label==1
        filename = sprintf('%s/%s%0.4d.tif',ImageFolder_label,'label',i); % Generate FILENAME of frame
        imwrite(uint16(im_label),filename,'tif'); % Write out tiff file
    end
    im_msr=measure(im_label,imdata,{'dimension','DimensionsEllipsoid','gravity','size','Mass','Mean','StdDev','MaxVal','MajorAxes'}, ...
        [],2,MinSpotSize*2,size(im,1)*size(im,2));
    clear dety2 detz2 eta_center MajorAxisAngle;
    for n=1:length(im_msr)   
        dety2(n)=ImageSize(1)-round(im_msr(n).Gravity(1));
        detz2(n)=ImageSize(2)-round(im_msr(n).Gravity(2));
        % angle with respect to vertical axis in coordinate system with origin sitting at the detector center
        % [0 2*pi]
        eta=atan2(-im_msr(n).Gravity(2)+ImageSize(2)/2,im_msr(n).Gravity(1)-ImageSize(1)/2); % [0, -pi] and [0, pi]
        if eta>0
            eta=pi/2-eta;
            if eta<0
                eta=eta+2*pi;
            end
        else
            eta=eta+3/2*pi;
        end
        eta=eta*180/pi;
        eta_center(n)=eta;
        MajorAxisAngle(n)=atand(im_msr(n).MajorAxes(2)/im_msr(n).MajorAxes(1));
        Spots{i}=[Spots{i};i n rotation_angle(i) eta dety2(n) detz2(n) im_msr(n).CartesianBox' ... % [1-8] Projtion no., labelID, rotation, eta, dety2, detz2, CartBox1, CartBox2 
            im_msr(n).Gravity' im_msr(n).Size' im_msr(n).MaxVal' im_msr(n).Mean' im_msr(n).StdDev' ... % [9 14] COM, size, IntMax, IntMean, IntSD
            im_msr(n).Mass' im_msr(n).DimensionsEllipsoid' ... % [15 17] IntegratedIntensity, EllipSoidaxis1, EllipSoidaxis2,
            tand(im_msr(n).MajorAxes(2)/im_msr(n).MajorAxes(1))]; % [18] major axis angle with respect to horizontal axis
    end
    varNames={'ProjNo','ObjectID','RotationAngle','eta','dety2','detz2','CartBox2', ...
        'COM','Size','IntMax','IntMean','IntSD','IntegratedIntensity','EllipSoidaxis', ...
        'MajorAxisAngle'};
    SpotsTable{i}=table(i*ones(length(im_msr),1),[1:length(im_msr)]',rotation_angle(i)*ones(length(im_msr),1), ...
        eta_center',dety2',detz2',im_msr.CartesianBox', ...
        im_msr.Gravity',im_msr.Size',im_msr.MaxVal',im_msr.Mean',im_msr.StdDev', ...
        im_msr.Mass',im_msr.DimensionsEllipsoid',MajorAxisAngle','VariableNames',varNames);
%     Sum - sum of object intensity (=mass) *
%     Mass - mass of object (=sum of object intensity) *
%     Mean - mean object intensity *
%     StdDev - standard deviation of object intensity *
%     MaxVal - maximum object intensity *
%     MajorAxes - principal axes of an object

    sprintf('get spot info for projection no. %d / %d',i,ImageNr)
end


%%%%%%%%%%%%%% add on Oct 5, 2021
for i=1:length(proj)
    proj_temp{i}=uint16(proj{i});
end
clear proj;
proj=proj_temp;
for i=1:length(proj_bin)
    proj_temp{i}=double(proj_bin{i});
end
clear proj_bin;
proj_bin=proj_temp;
clear proj_temp;
save(OutputFileName,'proj','proj_bin','Spots','SpotsTable','rot_start','rot_end','rot_step','ImageNr','-v7.3');



