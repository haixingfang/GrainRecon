% makeimage script produces tiff diffraction images using the reflection information in A
% Haixing Fang March 2019
% used for produces tiff diffraction images for polychromatic X-ray diffraction
% modified on Sep 25, 2019: use the searching method to find unique
% reflections to avoid connecting issues caused by using 'measure'
% employ Gaussian point spread
% last updated on March 30, 2020
% require DIPimage toolbox
% modified on Sep 7, 2021
peakshape=1;
peakfwhm=2;
pixellimit = ceil(peakfwhm);
label_unique_flag=1;

nr = size(A,1);
ImStack_frame_cl=newcolorim([detysize detzsize],'RGB','uint8');
ImStack_frame_cl_mask=newcolorim([detysize detzsize],'RGB','uint8');
ImStack_frame_cl_mask_rest=~newcolorim([detysize detzsize],'RGB','uint8');
nrefl=0;
Arefl=zeros(size(A));
frame = zeros(detzsize,detysize);
for jj=1:nr
    int=A(jj,21);
    % Changed such that dety,detz has 0,0 in the center of the
    % lower right corner instead of 0,0 at the border.
%     dety=round(A(jj,17))+1;
%     detz=round(A(jj,18))+1;
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

frame=ceil(frame);
frame = frame + bgint*ones(detzsize,detysize); % Add constant background counts
frame_image=uint16(frame);
frame_bin = frame_image>thres1;
if thres1~=(bgint+sqrt(bgint))
    frame_image=double(frame_image).*double(frame_bin)+thres1.*double(~frame_bin);
    frame_image=uint16(frame_image);
end
% add beam stop
if ~isempty(BeamStopY) && ~isempty(BeamStopZ)
    frame_image_BeamStop=frame_image;
    frame_image_BeamStop(BeamStopZ(1):BeamStopZ(2),BeamStopY(1):BeamStopY(2))=0;      
    frame_bin_BeamStop=frame_bin;
    frame_bin_BeamStop(BeamStopZ(1):BeamStopZ(2),BeamStopY(1):BeamStopY(2))=0; 
end
label_BeamStop=1;
if label_BeamStop~=1
    frame_label = label(frame_bin,1,0,detzsize*detysize);
    frame_label_im=uint16(frame_label);
    frame_measure = measure(frame_label,frame_image,{'dimension','DimensionsCube','gravity','size','feret'});  
else
    frame_label = label(frame_bin_BeamStop,1,0,detzsize*detysize);
    frame_label_im=uint16(frame_label);
    frame_measure = measure(frame_label,frame_image,{'dimension','DimensionsCube','gravity','size','feret'}); 
end

clear GrainIndex;
clear label_str;
label_pos=[];
Gr_unique=unique(A(:,2),'rows');
n=0;     % spot number
n_pos=0; % number for label string
for m=1:length(Gr_unique)
    hklIndex=A(find(A(:,2)==Gr_unique(m)),4:6);
    [hklIndex,ia,ic] = unique(hklIndex,'rows');
    hklIndex_unique=[];
    for n1=1:length(hklIndex(:,1))
        hkl_rest=setdiff(hklIndex,hklIndex(n1,:),'rows');
        if ~ismember(hklIndex(n1,:)/2,hkl_rest,'rows') && ~ismember(hklIndex(n1,:)/3,hkl_rest,'rows') ...
            && ~ismember(hklIndex(n1,:)/4,hkl_rest,'rows') && ~ismember(hklIndex(n1,:)/5,hkl_rest,'rows')
            hklIndex_unique=[hklIndex_unique;hklIndex(n1,:)]; % get rid of multiple reflections such as 111 222
        end
    end
    hklIndex_unique(:,4)=hklIndex_unique(:,1).^2+hklIndex_unique(:,2).^2+hklIndex_unique(:,3).^2;
    hklIndex_unique = sortrows(hklIndex_unique,4);
    hklIndex_unique=hklIndex_unique(:,1:3);
%         m
%         length(hklIndex_unique(:,1))
    for nn=1:length(hklIndex_unique(:,1))
        clear A_filted;
        A_filted=A(find(A(:,2)==Gr_unique(m) & ((A(:,4)==hklIndex_unique(nn,1) & ...
            A(:,5)==hklIndex_unique(nn,2) & A(:,6)==hklIndex_unique(nn,3)) | ...
            (A(:,4)==hklIndex_unique(nn,1)*2 & ...
            A(:,5)==hklIndex_unique(nn,2)*2 & A(:,6)==hklIndex_unique(nn,3)*2) | ...
            (A(:,4)==hklIndex_unique(nn,1)*3 & ...
            A(:,5)==hklIndex_unique(nn,2)*3 & A(:,6)==hklIndex_unique(nn,3)*3))),:);
        if ~isempty(A_filted)
        if label_BeamStop==1 && ~(round(min(A_filted(:,18)))>=BeamStopZ(1) && round(max(A_filted(:,18)))<=BeamStopZ(2) ...
                && round(min(A_filted(:,17)))>=BeamStopY(1) && round(max(A_filted(:,17)))<=BeamStopY(2))
            n=n+1;
            A_rest=setdiff(A,A_filted,'rows');
            if exist('ismembertol','builtin')==0
                LIA=ismember(A_filted(:,17:18),A_rest(:,17:18),'Rows');
            else
                LIA=ismembertol(A_filted(:,17:18),A_rest(:,17:18),2,'DataScale',1,'ByRows',true); % tolerance is 2 pixels
            end
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
%                 frame_select=imcrop(frame_image,[CropROI(1) CropROI(2) CropROI(3)-1 CropROI(4)-1]); % [xmin ymin width height]
            frame_select=imcrop(frame_image_BeamStop,[CropROI(1) CropROI(2) CropROI(3)-1 CropROI(4)-1]); % [xmin ymin width height]                
            if sum(sum(frame_select==0))>0 || nansum(nansum(frame_select==mean(nanmean(frame_select))))/numel(frame_select)
                frame_select_bg=1.1*thres1;
                frame_select_bin = frame_select>frame_select_bg; % in touch with the beam stop
            else
% %                     [frame_select_bin,frame_select_bg] = threshold(frame_select,'otsu',Inf);
%                     [frame_select_bin,frame_select_bg] = threshold(frame_select,'background',Inf);
					if sum(LIA)==0 % only one spot in the ROI
%                         if max(max(frame_select))<=thres1*1.5
                        [frame_select_bin,frame_select_bg1] = threshold(frame_select,'otsu',Inf);
%                         else
                        [frame_select_bin,frame_select_bg2] = threshold(frame_select,'background',Inf);
%                         end
						frame_select_bg=1/2*frame_select_bg1+1/2*frame_select_bg2;
					else
						[frame_select_bin,frame_select_bg] = threshold(frame_select,'otsu',Inf);
					end
                frame_select_bin=frame_select>frame_select_bg;
                if frame_select_bg<thres1
                    frame_select_bin=frame_select>thres1;
                end
                frame_select_bin=double(frame_select_bin);
            end
            frame_select_bin=frame_select_bin>0;
%                 dipshow(frame_select_bin);
            if sum(LIA)==0 % no overlap with other spots
%                    frame_select=imcrop(frame_image,[round(min(A_filted(:,17))) round(min(A_filted(:,18))) ...
%                        round(max(A_filted(:,17))-min(A_filted(:,17))) ...
%                        round(max(A_filted(:,18))-min(A_filted(:,18)))]); % [xmin ymin width height]
%                    frame_select_bin = frame_select>thres1;
%     %                frame_select_label = label(frame_select_bin,1,0,size(frame_select,1)*size(frame_select,2));
%     %                frame_select_measure = measure(frame_select_label,frame_select,{'gravity','size'});

               GrainIndex(n,1)=n; % ID of diffraction spot
%                GrainIndex(n,2)=sum(frame_select_measure.Size); % size of the diffraction spot [pixel*pixel]
               GrainIndex(n,2)=sum(sum(frame_select_bin));
               GrainIndex(n,3)=A_filted(1,2); % ID of grain
               GrainIndex(n,4:6)=hklIndex_unique(nn,:); % (h k l)
               GrainIndex(n,7)=sum(sum(double(frame_select_bin).*double(frame_select)))/GrainIndex(n,2);% average intensity of the diffraction spot
               GrainIndex(n,8)=A_filted(1,1); % ID of Reflection number
               GrainIndex(n,9)=0; % overlapped area fraction, flag 0-no, 1-yes
               GrainIndex(n,10)=sum(sum(double(frame_select_bin).*double(frame_select)))- ...
                   sum(sum(double(frame_select_bin).*max([thres1 frame_select_bg]))); % IntInt, modified on March 28, 2020
               AllSpotsNr(rot_number)=AllSpotsNr(rot_number)+1;
            else           % overlap with other spots
%                    frame_select=imcrop(frame_image,[round(min(A_filted(:,17))) round(min(A_filted(:,18))) ...
%                        round(max(A_filted(:,17))-min(A_filted(:,17))) ...
%                        round(max(A_filted(:,18))-min(A_filted(:,18)))]); % [xmin ymin width height]
%                    frame_select_bin = frame_select>thres1;
%     %                frame_select_label = label(frame_select_bin,1,0,size(frame_select,1)*size(frame_select,2));
%     %                frame_select_measure = measure(frame_select_label,frame_select,{'gravity','size'});

               GrainIndex(n,1)=n; % ID of diffraction spot
%                GrainIndex(n,2)=sum(frame_select_measure.Size); % size of the diffraction spot [pixel*pixel]
               GrainIndex(n,2)=sum(sum(frame_select_bin));
               GrainIndex(n,3)=A_filted(1,2); % ID of grain
               GrainIndex(n,4:6)=hklIndex_unique(nn,:); % (h k l)
               GrainIndex(n,7)=sum(sum(double(frame_select_bin).*double(frame_select)))/GrainIndex(n,2);% average intensity of the diffraction spot
               GrainIndex(n,8)=A_filted(1,1); % ID of Reflection number
               GrainIndex(n,9)=length(find(LIA==1))/length(LIA); % overlapped area fraction, flag 0-no, 1-yes
               GrainIndex(n,10)=sum(sum(double(frame_select_bin).*double(frame_select)))- ...
                   sum(sum(double(frame_select_bin).*max([thres1 frame_select_bg]))); % IntInt, modified on March 28, 2020
               AllSpotsNr(rot_number)=AllSpotsNr(rot_number)+1;
               SpotOverlapNr(rot_number)=SpotOverlapNr(rot_number)+1;
            end
            if GrainIndex(n,2)>0
                n_pos=n_pos+1;
                label_str{n_pos} = strcat([num2str(GrainIndex(n,4)) num2str(GrainIndex(n,5)) ...
                    num2str(GrainIndex(n,6))],strcat(', No.',num2str(GrainIndex(n,3))));
                label_pos = [label_pos;mean(A_filted(:,17)) mean(A_filted(:,18))];
%                 label_pos = [label_pos;detysize-mean(A_filted(:,17)) detzsize-mean(A_filted(:,18))];
            end
			hklno=find(hkl_square==(hklIndex_unique(nn,1)^2+hklIndex_unique(nn,2)^2+hklIndex_unique(nn,3)^2));
%                 frame_select_bin_edge=dgg(frame_select_bin,1);
%                 frame_select_bin_edge=frame_select_bin_edge~=0;
%                 frame_select_bin_edge=berosion(frame_select_bin_edge,3,-2,1);
			im1 = countneighbours(frame_select_bin,1,Inf,0);
			im2=im1>=8;
			frame_select_bin_edge=frame_select_bin-im2;
%                 frame_select_bin_edge=bclosing(frame_select_bin_edge,2,-1,0);
			frame_select_bin_edge=bdilation(frame_select_bin_edge,1,-1,0);
			frame_cl=colorspace(frame_select_bin_edge,'RGB');
			frame_cl_mask=~frame_cl;
			frame_cl=dip_array(frame_cl,'double');
			redChannel = frame_cl(:, :, 1)*hkl_color(hklno,1);
			greenChannel = frame_cl(:, :, 2)*hkl_color(hklno,2);
			blueChannel = frame_cl(:, :, 3)*hkl_color(hklno,3);
			frame_cl=joinchannels('RGB',redChannel, greenChannel, blueChannel);            
			ImStack_frame_cl(CropROI(1):CropROI(1)+CropROI(3)-1,CropROI(2):CropROI(2)+CropROI(4)-1)=frame_cl;
			ImStack_frame_cl_mask(CropROI(1):CropROI(1)+CropROI(3)-1,CropROI(2):CropROI(2)+CropROI(4)-1)=frame_cl_mask;
			frame_cl_mask_rest=newcolorim([length(frame_select_bin(1,:)) length(frame_select_bin(:,1))],'RGB','uint8');
			ImStack_frame_cl_mask_rest(CropROI(1):CropROI(1)+CropROI(3)-1,CropROI(2):CropROI(2)+CropROI(4)-1)=frame_cl_mask_rest;
        end
        end
    end
end
% unique label annotations
if exist('label_str','var')
    if label_unique_flag==1
        [label_str_unique,idx,idc]=unique(label_str);
        label_pos_unique=label_pos(idx,:);
    else
        label_str_unique=label_str;
        label_pos_unique=label_pos;
    end
end

% flipud(fliplr) is to rotate the image by 180 degrees in clockwise direction

%Write out tiff file for grey value image
filename = sprintf('%s/%s%0.4d_grey.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
% frame_image = flipud(fliplr(frame_image)); %flip frame to output images in correct direction
%     if rot<0
%         frame_image = fliplr(frame_image);
%     end
% imwrite(frame_image,filename,'tif'); % Write out tiff file

%Write out tiff file for binary image
filename = sprintf('%s/%s%0.4d_binary.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
% frame_bin = flipud(fliplr(frame_bin)); %flip frame to output images in correct direction
%     if rot<0    
%         frame_bin = fliplr(frame_bin);
%     end
% imwrite(frame_bin,filename,'tif'); % Write out tiff file

%Write out tiff file for label image
filename = sprintf('%s/%s%0.4d_label.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
% frame_label_im = flipud(fliplr(frame_label_im)); %flip frame to output images in correct direction
% imwrite(frame_label_im,filename,'tif'); % Write out tiff file

Im=double(ImStack_frame_cl);
% Im_fliplr = flip(Im,2);
% Im_flipudlr = flip(Im_fliplr,1);
Im_flipudlr=Im;
redChannel = Im_flipudlr(:, :, 1);
greenChannel = Im_flipudlr(:, :, 2);
blueChannel = Im_flipudlr(:, :, 3);
ImStack_frame_cl=joinchannels('RGB',redChannel, greenChannel, blueChannel);

Im=double(ImStack_frame_cl_mask);
% Im_fliplr = flip(Im,2);
% Im_flipudlr = flip(Im_fliplr,1);
Im_flipudlr=Im;
redChannel = Im_flipudlr(:, :, 1);
greenChannel = Im_flipudlr(:, :, 2);
blueChannel = Im_flipudlr(:, :, 3);
ImStack_frame_cl_mask=joinchannels('RGB',redChannel, greenChannel, blueChannel);

Im=double(ImStack_frame_cl_mask_rest);
% Im_fliplr = flip(Im,2);
% Im_flipudlr = flip(Im_fliplr,1);
Im_flipudlr=Im;
redChannel = Im_flipudlr(:, :, 1);
greenChannel = Im_flipudlr(:, :, 2);
blueChannel = Im_flipudlr(:, :, 3);
ImStack_frame_cl_mask_rest=joinchannels('RGB',redChannel, greenChannel, blueChannel);
ImStack_frame_cl_mask=ImStack_frame_cl_mask+ImStack_frame_cl_mask_rest;

%Write out tiff file for label image
if exist('label_pos_unique','var') && exist('insertText.m','file')~=0
    if ~isempty(label_pos_unique)
        filename = sprintf('%s/%s%0.4d_label_annot.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
    %     frame_label_im = flipud(fliplr(image_label_annot)); %flip frame to output images in correct direction
        frame_label_annot = insertText(frame_label_im, label_pos_unique, label_str_unique, 'FontSize', 30,'BoxOpacity', 0.4);
        imwrite(frame_label_annot,filename,'tif'); % Write out tiff file
    end
end

%Write out tiff file for grey value image with beam stop
filename = sprintf('%s/%s%0.4d_beamstop.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
% frame_image_BeamStop = flipud(fliplr(frame_image_BeamStop));
imwrite(frame_image_BeamStop,filename,'tif'); % Write out tiff file

%Write out tiff file for binary image with beam stop
filename = sprintf('%s/%s%0.4d_beamstop_bin.tif',direc,prefix,rot_number-1); % Generate FILENAME of frame
% frame_bin_BeamStop = flipud(fliplr(frame_bin_BeamStop));
imwrite(frame_bin_BeamStop,filename,'tif'); % Write out tiff file

% record the effective index of diffraction spots
if exist('GrainIndex','var')
    hklIndex=GrainIndex(:,3:6);
    [hklIndex_unique,ia,ic] = unique(hklIndex,'rows');
    GrainIndex_unique=GrainIndex(find(GrainIndex(:,2)>0),:);
end


