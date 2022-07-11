% spot segmentation by combining Laplacian of Gaussian and single
% thresholding

function [proj proj_bin Spots SpotsTable rot_start rot_end rot_step ImageNr]=SpotsFromBinaryImage(ImageFolder,name1,name3, ...
    Nrolling,Background,sigma,percent,MinSpotSize,CombineFlag,write_im,parallelflag)

% for testing
% clear all;
% ImageFolder='M:\MetaData\Projections_Dss14_Dsd14\Projections_16';
% name1='proj';
% name3='.tiff';
% Nrolling=11;
% BackGround=5; % [0 6] 
% sigma=6; % [1-25]
% percent=0.10; % [0 0.25] minimum intensity fraction of a blob to be considered as a spot
% MinSpotSize=10; % 2-10 [pixels]
% CombineFlag=1;
% write_im=1; % save segmented image
% parallelflag=1; % parallel computing, only suitable to images with many spots

if write_im==1 % write images
    if ~exist('.\SpotIm','dir')
        mkdir('.\SpotIm');
        OutputFolder='.\SpotIm';
    else
        OutputFolder='.\SpotIm';
    end
end
if exist(ImageFolder, 'dir')~=7
  Message = sprintf('Error: The following folder does not exist:\n%s', ImageFolder)
%  uiwait(warndlg(Message));
  return;
end
ImageInfo=dir([ImageFolder '/*' name3]);
ImageNr = numel(ImageInfo);
rot_start=-180;
rot_end=180;
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
    if rot_number==1
        detysize=imsize(proj{1},1);
        detzsize=imsize(proj{1},2);
        if detysize>2000 % no binning
            BeamStopY=[560 1455];
            BeamStopZ=[545 1460];
        elseif detysize>1000 && detysize<1100 % binning 2
            BeamStopY=[294 723];
            BeamStopZ=[278 725];
        else
            BeamStopY=[560 1455];
            BeamStopZ=[545 1460];
            sprintf('Warning: check the position of the beam stop first!!!')
        end
    end
    sprintf('Reading projection no. %d / %d',rot_number,ImageNr)
    rot_number=rot_number+1;
end
% Nrolling=11; % recommended values are 10-15 in GrainMapper
% Nrolling=round(ImageNr/3);
rev_flag=1;
if rev_flag==0
    rolling=Rolling_Median(proj,Nrolling);
    for i=1:length(proj)
        j=fix(i/Nrolling)+1;
        if j>length(rolling)
            j=length(rolling);
        end
        proj_roll{i}=proj{i}-rolling{j}; % subtract the rolling median image, regarded as subtract the Bg
        sprintf('Rolling median fitering process of projection no. %d / %d',i,ImageNr)
    end
else
    rolling=Rolling_Median_rev(proj,Nrolling);
    for i=1:length(proj)
        proj_roll{i}=proj{i}-rolling{i}; % subtract the rolling median image, regarded as subtract the Bg
        sprintf('Rolling median fitering process of projection no. %d / %d',i,ImageNr)
    end
end

% new algorithm set up on Sep 4, 2020
%  Jonathan F. Lind,
% "In-situ High-Energy diffraction Microscopy Study of Zirconium Under Uni-axial Tensile Deformation"
% 2013, PhD thesis, page 20

if CombineFlag==1
    ThresholdValue=125; % binary segmentation
end
for i=1:length(proj_roll)
    Spots{i}=[];
    im=proj_roll{i};
    if i==1
       ImageSize=size(proj_roll{i});
    end
    im_spot_bin=LoG_segmentation(im,BackGround,sigma,percent,MinSpotSize,parallelflag,BeamStopY,BeamStopZ);
%     im_spot_bin(BeamStopY(1):BeamStopY(2),BeamStopZ(1):BeamStopZ(2))=0;
    im_spot_bin=dip_image(im_spot_bin,'bin');
    if CombineFlag==1
        im_bin=im>ThresholdValue;
        im_spot_bin=im_spot_bin+im_bin;
        im_spot_bin=im_spot_bin>0;
    end
    proj_bin{i}=im_spot_bin;
    proj_spot{i}=im_spot_bin*double(im);
    imdata=proj{i};
    im_label=label(im_spot_bin,2,MinSpotSize*2,size(im,1)*size(im,2));
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

    if write_im==1 % write images
        filename = sprintf('%s/%s%0.4d.tif',OutputFolder,'proj',i); % Generate FILENAME of frame
        frame=uint8(im_spot_bin*255);
        imwrite(frame,filename,'tif');
    end
    sprintf('LoG filtering and spot segmentation for projection no. %d / %d',i,ImageNr)
end


