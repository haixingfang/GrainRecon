function [GrainIndex,GrainIndex_unique]=make_projection_simu(A,SubGrain, ...
                Lsam2sou,Lsam2det,pixelysize,pixelzsize,detysize,detzsize, ...
                direc,prefix,rot_number,BeamStopY,BeamStopZ)
                
% make_projection script produces tiff diffraction images using the reflection information in A
% used for produces tiff diffraction images for polychromatic X-ray diffraction
% modified on Sep 25, 2019: use the searching method to find unique
% reflections to avoid connecting issues caused by using 'measure'
% employ Gaussian point spread
% last updated on March 30, 2020
% require DIPimage toolbox
% modified on Sep 7, 2021
% does not require DIPimage toolbox
    peakshape=1;
    peakfwhm=2;
    pixellimit = ceil(peakfwhm);
    label_unique_flag=1;

    nr = size(A,1);
    ImStack_frame_cl=zeros([detzsize detysize 3],'uint8');
    ImStack_frame_cl_mask=zeros([detzsize detysize 3],'uint8');
    ImStack_frame_cl_mask_rest=uint8(~ImStack_frame_cl_mask);

    nrefl=0;
    Arefl=zeros(size(A));
    frame = zeros(detzsize,detysize);
    for jj=1:nr
        int=A(jj,21);
        % Changed such that dety,detz has 0,0 in the center of the
        % lower right corner instead of 0,0 at the border.
%         dety=round(A(jj,17))+1;
%         detz=round(A(jj,18))+1;
        % After May 12, 2022
        dety=round(A(jj,17));
        detz=round(A(jj,18));

        % Do not consider reflections with a CMS further away from the
        % detector than 5*fwhm of the spot
        if (-5*pixellimit >= dety) || (dety >= detysize+5*pixellimit) || (-5*pixellimit >= detz) || (detz >= detzsize+5*pixellimit)
            %disp(['reflection outside detector; y,z: ',num2str(dety),', ',num2str(detz)])
        else
            nrefl=nrefl+1;
            Arefl(nrefl,:) = A(jj,:);
            if peakshape == 0                 % spike peak
                pixelnr_fit(1)=0.0943; % 14-14, calibrated
                pixelnr_fit(2)=1.1679; % 14-14
                if SubGrain{A(jj,2)}(A(jj,23),6)>1 && SubGrain{A(jj,2)}(A(jj,23),6)<100 % use um
    %                         pixelnr=fix(SubGrain{A(jj,2)}(A(jj,23),6)./(2*1000*mean([pixelysize pixelzsize]))); % [2 1.75 1.5 1.25]
                    pixelnr=SubGrain{A(jj,2)}(A(jj,23),6)*pixelnr_fit(1)+pixelnr_fit(2); % 14-14
                elseif SubGrain{A(jj,2)}(A(jj,23),6)==Inf
                    pixelnr=5;
                else
    %                         pixelnr=fix(SubGrain{A(jj,2)}(A(jj,23),6)./(2*mean([pixelysize pixelzsize]))); % use mm
                    pixelnr=SubGrain{A(jj,2)}(A(jj,23),6)*1000*pixelnr_fit(1)+pixelnr_fit(2); % 14-14
                end
                pixelnr=pixelnr*(Lsam2det/Lsam2sou);
    %                     pixelnr=3;

    %                     % randomize the position to avoid 'jaggie' artifacts
                dety_delta=randi([dety-round(pixelnr/2),dety+round(pixelnr/2)],1,1)-dety;
                dety = dety+dety_delta;
                rand_direction=rand(1);
                if rand_direction>0.5
                    detz = detz+dety_delta*tand(26.6); % rotated grid method
                else
                    detz = detz-dety_delta*tand(26.6); % rotated grid method
                end
                for k1=round(dety-pixelnr):round(dety+pixelnr)
                    for k2=round(detz-pixelnr):round(detz+pixelnr)
                        if (0 < k1) && (k1 <= detysize) && (0 < k2) && (k2 <= detzsize)
                            frame(k2,k1)=frame(k2,k1)+int;
                        end
                    end
                end
            else             % Gaussian type peak
                % find factor of peak on this frame
                factor=1/(1/(2*pi*4*4)*exp(-1))*0.1; % a factor for intensity for adjusting appearing range, no physical meaning
                if SubGrain{A(jj,2)}(A(jj,23),6)>1 && SubGrain{A(jj,2)}(A(jj,23),6)<100 % use um
                    pixelnr=(SubGrain{A(jj,2)}(A(jj,23),6)./(2*1000*mean([pixelysize pixelzsize]))); % [2 1.75 1.5 1.25]
                elseif SubGrain{A(jj,2)}(A(jj,23),6)==Inf
                    pixelnr=5;
                else
                    pixelnr=(SubGrain{A(jj,2)}(A(jj,23),6)./(2*mean([pixelysize pixelzsize]))); % use mm
                end
                pixelnr=pixelnr*(Lsam2det/Lsam2sou);
                if pixelnr<=0
                    pixelnr=1;
                end
                peakfwhm=sqrt(2)*pixelnr;
                pixelnr=9/20*2.355*peakfwhm;% FWHM/FWTM = 5/9
                gaussian1=fspecial('Gaussian',round(2*pixelnr)+1,peakfwhm);
                PSF2=fspecial('motion',pixelnr*2,90-abs(A(jj,16)-90)); % anisotropic filter
                PSF2=conv2(PSF2,gaussian1); % convolution of 'motion' and Gaussian filters
                for k1=1:size(PSF2,2)
                    for k2=1:size(PSF2,1)
                        if (0 < round(k1+dety-size(PSF2,2)/2)) && (round(k1+dety-size(PSF2,2)/2) <= detysize) ...
                                && (0 < round(k2+detz-size(PSF2,1)/2)) && (round(k2+detz-size(PSF2,1)/2) <= detzsize)
                            frame(round(k2+detz-size(PSF2,1)/2),round(k1+dety-size(PSF2,2)/2))= ...
                                frame(round(k2+detz-size(PSF2,1)/2),round(k1+dety-size(PSF2,2)/2))+int*factor*PSF2(k2,k1);
                        end
                    end
                end
            end    
        end
    end

%     frame_image=uint16(ceil(frame));
    frame_image=frame;clear frame;
    ThresholdValue=0.1; % 0.2? if it is 0, it is possible to segment spots as one sigle spot
    frame_bin = frame_image>ThresholdValue;
%     [SpotBoundary,frame_label,~,~] = bwboundaries(frame_bin,'noholes');
    frame_label=bwlabeln(frame_bin,18);
    % frame_measure=regionprops(frame_label,'all');
    
    %{
    % add on Nov 8
    percent=0.15;
    im_spot=zeros(size(frame_image));
    for j=1:max(max(frame_label))
        ConnectedID=frame_label==j;
        ConnectedComp=ConnectedID.*frame_image;
        thres_int=percent*max(max(ConnectedComp));
        AcceptComp=ConnectedComp>thres_int;
        im_spot=im_spot+AcceptComp;
        j;
    end
    im_spot_bin=im_spot>0;
    frame_label=bwlabeln(im_spot_bin,18);
    %}
    
    frame_measure=regionprops(frame_label,frame_image,'Area','Centroid','BoundingBox', ...
                'EquivDiameter','MeanIntensity','PixelValues','WeightedCentroid');
    frame_measure_list=zeros(length(frame_measure),8);
    for kk=1:length(frame_measure)
        frame_measure_list(kk,:)=[frame_measure(kk).Area frame_measure(kk).WeightedCentroid frame_measure(kk).BoundingBox ...
                                   frame_measure(kk).EquivDiameter];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% summary
    clear GrainIndex;
    clear label_str;
    A_unique=unique(A(:,[2 4:6]),'rows');
    n=0;
    n_pos=0;
    se = strel('sphere',1);
    label_pos=[];
    SpotsPair=[];
    if ~isempty(frame_measure_list)
    for nn=1:length(A_unique(:,1))
        A_filted=A(A(:,2)==A_unique(nn,1) & ((A(:,4)==A_unique(nn,2) & ...
            A(:,5)==A_unique(nn,3) & A(:,6)==A_unique(nn,4)) | ...
            (A(:,4)==A_unique(nn,2)*2 & ...
            A(:,5)==A_unique(nn,3)*2 & A(:,6)==A_unique(nn,4)*2) | ...
            (A(:,4)==A_unique(nn,2)*3 & ...
            A(:,5)==A_unique(nn,3)*3 & A(:,6)==A_unique(nn,4)*3)),:);       
        dis_all=sqrt((mean(A_filted(:,17))-frame_measure_list(:,2)).^2+ ...
            (mean(A_filted(:,18))-frame_measure_list(:,3)).^2);
        [dis_all_min,pair_spot_ID]=min(dis_all);
		if dis_all_min<15 % otherwise it matches with a wrong segmented forward spot, which may not be segmented
        n=n+1;
        GrainIndex(n,1)=pair_spot_ID; % ID of the spot
        GrainIndex(n,2)=frame_measure_list(pair_spot_ID,1); % size of the spot
        GrainIndex(n,3)=A_filted(1,2); % ID of the grain
        GrainIndex(n,4:6)=A_filted(1,4:6); % hkl
        GrainIndex(n,7)=frame_measure(pair_spot_ID).MeanIntensity; % mean intensity
%         GrainIndex(n,8)=A_filted(1,1); % ID of reflection number
        GrainIndex(n,8)=A_filted(1,22); % spot energy
        GrainIndex(n,9)=0; % overlapped area fraction, flag 0-no, 1-yes
        GrainIndex(n,10)=sum(frame_measure(pair_spot_ID).PixelValues); % integrated intensity
        GrainIndex(n,11)=A_filted(1,14); % rot

%         if GrainIndex(n,2)>0
            n_pos=n_pos+1;
            label_str{n_pos} = strcat([num2str(GrainIndex(n,4)) num2str(GrainIndex(n,5)) ...
                num2str(GrainIndex(n,6))],strcat(', No.',num2str(GrainIndex(n,3))));
            label_pos = [label_pos;mean(A_filted(:,17)) mean(A_filted(:,18))];
%         end

        %{
        % pairing the spots
        spotpair_exp_ID=find(sqrt((Spots{rot_number}(:,9)-frame_measure_list(pair_spot_ID,2)).^2+ ...
            (Spots{rot_number}(:,10)-frame_measure_list(pair_spot_ID,3)).^2)<15 & ...
            (Spots{rot_number}(:,11)-frame_measure_list(pair_spot_ID,1))./frame_measure_list(pair_spot_ID,1)<5 & ...
            (Spots{rot_number}(:,11)-frame_measure_list(pair_spot_ID,1))./frame_measure_list(pair_spot_ID,1)>-0.5);
        if ~isempty(spotpair_exp_ID) && A_filted(1,4)^2+A_filted(1,5)^2+A_filted(1,6)^2<=hkl_square(3) %...
            % && A_filted(1,22)>=7 && A_filted(1,22)<=25 % spot energy should be in the range [7 30]
            if length(spotpair_exp_ID)>1
                [~,ind]=min(sqrt((Spots{rot_number}(spotpair_exp_ID,9)-frame_measure_list(pair_spot_ID,2)).^2+ ...
                    (Spots{rot_number}(spotpair_exp_ID,10)-frame_measure_list(pair_spot_ID,3)).^2));
                spotpair_exp_ID=spotpair_exp_ID(ind);
            end
            % grainID,pos,hkl,euler_angles,energy_hkl,rot_angle,dety,detz,spotpair_exp_ID,dety_exp,detz_exp,disy,disz,dis
            SpotsPair=[SpotsPair;A_filted(1,2) mean(SubGrain{A_filted(1,2)}(:,2:4),1) A_filted(1,4:6) ...
                    A_filted(1,8:10) A_filted(1,22) rot ...
                    frame_measure_list(pair_spot_ID,2:3) spotpair_exp_ID Spots{rot_number}(spotpair_exp_ID,9:10) ...
                    frame_measure_list(pair_spot_ID,2)-Spots{rot_number}(spotpair_exp_ID,9) ...
                    frame_measure_list(pair_spot_ID,3)-Spots{rot_number}(spotpair_exp_ID,10) ...
                    sqrt((Spots{rot_number}(spotpair_exp_ID,9)-frame_measure_list(pair_spot_ID,2)).^2+ ...
                        (Spots{rot_number}(spotpair_exp_ID,10)-frame_measure_list(pair_spot_ID,3)).^2) ...
                        frame_measure_list(pair_spot_ID,1) Spots{rot_number}(spotpair_exp_ID,11)];
        end
        
        hklno=find(hkl_square==(GrainIndex(n,4)^2+GrainIndex(n,5)^2+GrainIndex(n,6)^2));
        CropROI(1)=round(min(A_filted(:,17)))-20;
        CropROI(2)=round(min(A_filted(:,18)))-20;
        CropROI(3)=round(max(A_filted(:,17))-min(A_filted(:,17)))+40;
        CropROI(4)=round(max(A_filted(:,18))-min(A_filted(:,18)))+40;
        if CropROI(1)<1
            CropROI(1)=1;
        end
        if CropROI(2)<1
            CropROI(2)=1;
        end
        % [xmin ymin width height]
        if CropROI(1)+CropROI(3)>=detysize
            CropROI(3)=detysize-CropROI(1);
        end
        if CropROI(2)+CropROI(4)>=detzsize
            CropROI(4)=detzsize-CropROI(2);
        end
        frame_select_bin=imcrop(frame_bin,[CropROI(1) CropROI(2) CropROI(3)-1 CropROI(4)-1]); % [xmin ymin width height]
        if sum(sum(frame_select_bin))>0
            SpotBoundary = bwboundaries(frame_select_bin,'noholes');
            frame_select_bin_edge=zeros(size(frame_select_bin),'uint8');
            ind=sub2ind(size(frame_select_bin),SpotBoundary{1}(:,1),SpotBoundary{1}(:,2));
            frame_select_bin_edge(ind)=1;
            frame_select_bin_edge=imdilate(frame_select_bin_edge,se);
            frame_cl = cat(3, frame_select_bin_edge, frame_select_bin_edge, frame_select_bin_edge);

            redChannel = frame_cl(:, :, 1)*hkl_color(hklno,1);
            greenChannel = frame_cl(:, :, 2)*hkl_color(hklno,2);
            blueChannel = frame_cl(:, :, 3)*hkl_color(hklno,3);
            frame_cl=cat(3,redChannel,greenChannel,blueChannel);

            frame_select_bin_mask=uint8(~frame_select_bin_edge);
            frame_cl_mask = cat(3, frame_select_bin_mask, frame_select_bin_mask, frame_select_bin_mask);

            ImStack_frame_cl(CropROI(2):CropROI(2)+CropROI(4)-1,CropROI(1):CropROI(1)+CropROI(3)-1,:)=frame_cl;
            ImStack_frame_cl_mask(CropROI(2):CropROI(2)+CropROI(4)-1,CropROI(1):CropROI(1)+CropROI(3)-1,:)=frame_cl_mask;
            frame_cl_mask_rest=zeros([size(frame_select_bin) 3],'uint8');
            ImStack_frame_cl_mask_rest(CropROI(2):CropROI(2)+CropROI(4)-1,CropROI(1):CropROI(1)+CropROI(3)-1,:)=frame_cl_mask_rest;
        end
        %}
		end
    end
    end
    
    % unique label annotations
    if exist('label_str','var')
        if label_unique_flag==1
            [label_str_unique,idx,~]=unique(label_str);
            label_pos_unique=label_pos(idx,:);
        else
            label_str_unique=label_str;
            label_pos_unique=label_pos;
        end
    end
    
    %Write out tiff file for grey value image
    filename = sprintf('%s/%s%0.4d_grey.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
%     imwrite(frame_image,filename,'tif'); % Write out tiff file
    imwrite(uint16(double(frame_image)*65536/max(frame_image(:))),filename,'tif'); % Write out tiff file
    
    %Write out tiff file for binary image
    filename = sprintf('%s/%s%0.4d_binary.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
    imwrite(frame_bin,filename,'tif'); % Write out tiff file

    %Write out tiff file for label image
    filename = sprintf('%s/%s%0.4d_label.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
    % frame_label_im = flipud(fliplr(frame_label_im)); %flip frame to output images in correct direction
%     imwrite(frame_label_im,filename,'tif'); % Write out tiff file

%     ImStack_frame_cl_mask=ImStack_frame_cl_mask+ImStack_frame_cl_mask_rest;

    %Write out tiff file for label image
    if exist('label_pos_unique','var') && exist('insertText.m','file')~=0
        if ~isempty(label_pos_unique)
            filename = sprintf('%s/%s%0.4d_label_annot.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
            frame_label_annot = insertText(frame_label, label_pos_unique, label_str_unique, 'FontSize', 30,'BoxOpacity', 0.4);
            imwrite(frame_label_annot,filename,'tif'); % Write out tiff file
        end
    end

    % add on Nov 24, 2021
%     max(frame_image(:))
    bgint=926;
    frame_image_BG = double(frame_image)*(65536-bgint)/max(frame_image(:))+bgint*ones(detzsize,detysize); % Add constant background counts
    frame_image_BG=uint16(frame_image_BG);
%     max(frame_image_BG(:))
    % add beam stop
    if ~isempty(BeamStopY) && ~isempty(BeamStopZ)
        frame_image_BG(BeamStopZ(1):BeamStopZ(2),BeamStopY(1):BeamStopY(2))=0;       
    end
    %Write out tiff file for grey value image with a constant background
    filename = sprintf('%s/%s%0.4d_greyBG.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
    imwrite(frame_image_BG,filename,'tif'); % Write out tiff file
    
%     % Write out tiff file for overlaid image
%     imdata_cl=cat(3,uint8(imdata),uint8(imdata),uint8(imdata));
%     ImOverlay=imdata_cl.*ImStack_frame_cl_mask+ImStack_frame_cl;
%     filename = sprintf('%s/%s%0.4d_overlay.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
%     imwrite(ImOverlay,filename,'tif'); % Write out tiff file
    
    %{
    % write out spot centroid tracking tiff
    if ~isempty(SpotsPair)
        imdata_cl_track=imdata_cl;
        for i=1:length(SpotsPair(:,1))
            if (round(SpotsPair(i,14))>2 && round(SpotsPair(i,14))<detzsize-2 && ...
                round(SpotsPair(i,13))>2 && round(SpotsPair(i,13))<detysize-2 && ...
                round(SpotsPair(i,17))>2 && round(SpotsPair(i,17))<detzsize-2 && ...
                round(SpotsPair(i,16))>2 && round(SpotsPair(i,16))<detysize-2)
                imdata_cl_track(round(SpotsPair(i,14))-2:round(SpotsPair(i,14))+2, ...
                    round(SpotsPair(i,13))-2:round(SpotsPair(i,13))+2,1)=255;
                imdata_cl_track(round(SpotsPair(i,14))-2:round(SpotsPair(i,14))+2, ...
                    round(SpotsPair(i,13))-2:round(SpotsPair(i,13))+2,2)=0;
                imdata_cl_track(round(SpotsPair(i,14))-2:round(SpotsPair(i,14))+2, ...
                    round(SpotsPair(i,13))-2:round(SpotsPair(i,13))+2,3)=0;   % simulated spots

                imdata_cl_track(round(SpotsPair(i,17))-2:round(SpotsPair(i,17))+2, ...
                    round(SpotsPair(i,16))-2:round(SpotsPair(i,16))+2,1)=0;
                imdata_cl_track(round(SpotsPair(i,17))-2:round(SpotsPair(i,17))+2, ...
                    round(SpotsPair(i,16))-2:round(SpotsPair(i,16))+2,2)=0;
                imdata_cl_track(round(SpotsPair(i,17))-2:round(SpotsPair(i,17))+2, ...
                    round(SpotsPair(i,16))-2:round(SpotsPair(i,16))+2,3)=255; % exp spots
            end
        end
%         imwrite(imdata_cl_track,sprintf('%s/%s%0.4d_track.tif',direc,prefix,rot_number-1),'tif'); % Write out tiff file
    end
    %}
    
    % record the effective index of diffraction spots
    if exist('GrainIndex','var')
        hklIndex=GrainIndex(:,3:6);
%         [hklIndex_unique,ia,ic] = unique(hklIndex,'rows');
        GrainIndex_unique=GrainIndex(GrainIndex(:,2)>0,:);
	else
        GrainIndex=[];
        GrainIndex_unique=[];
    end
% end
% % visualization with DIPimage
% im=ImStack_frame_cl;
% im=ImStack_frame_cl_mask_rest;
% im1=joinchannels('RGB',im(:,:,1), im(:,:,2), im(:,:,3)); 
% dipshow(im1)
