function [Im_GB_rec_cl]=slice_show(DS,SliceNumber,compare_XY,compare_XZ,compare_YZ)

if compare_XY==1
    slice_rec=DS.GIDvol(:,:,SliceNumber);
    Im_GB_rec=(reshape(slice_rec,[DS.Dimension(1) DS.Dimension(2)]));
    Im_GB_rec=imrotate(Im_GB_rec,-90);
    Im_GB_rec=(flipud(fliplr(Im_GB_rec)));
    Dimension=[DS.Dimension(1) DS.Dimension(2)];
    redChannel = DS.IPF001(1,:,:,SliceNumber);
    redChannel=(reshape(redChannel,[DS.Dimension(1) DS.Dimension(2)])).*255;
    redChannel=imrotate(redChannel,-90);
    redChannel=(flipud(fliplr(redChannel)));
    greenChannel = DS.IPF001(2, :, :, SliceNumber);
    greenChannel=(reshape(greenChannel,[DS.Dimension(1) DS.Dimension(2)])).*255;
    greenChannel=imrotate(greenChannel,-90);
    greenChannel=(flipud(fliplr(greenChannel)));
    blueChannel = DS.IPF001(3, :, :, SliceNumber);
    blueChannel=(reshape(blueChannel,[DS.Dimension(1) DS.Dimension(2)])).*255;
    blueChannel=imrotate(blueChannel,-90);
    blueChannel=(flipud(fliplr(blueChannel)));
    Im_GB_rec_cl=joinchannels('RGB',redChannel, greenChannel, blueChannel);
    % add on Oct 1, 2020
    Im_GB_rec_bin=Im_GB_rec>0;
    Im_mask=~Im_GB_rec_bin;
    Im_mask=bwareaopen(Im_mask,round(size(Im_mask,1)*size(Im_mask,2)*0.01));
    Im_mask=Im_mask*255; % white background
    redChannel_mask=Im_mask;
    greenChannel_mask=Im_mask;
    blueChannel_mask=Im_mask;
    Im_mask_cl=joinchannels('RGB',redChannel_mask, greenChannel_mask, blueChannel_mask);
    Im_GB_rec_cl=Im_GB_rec_cl+Im_mask_cl;
elseif compare_XZ==1
    slice_rec=DS.GIDvol(:,SliceNumber,:);
    Im_GB_rec=(reshape(slice_rec,[DS.Dimension(1) DS.Dimension(3)]));
    Im_GB_rec=imrotate(Im_GB_rec,-90);
    Im_GB_rec=(flipud(fliplr(Im_GB_rec)));
    Dimension=[DS.Dimension(1) DS.Dimension(3)];
    redChannel = DS.IPF001(1, :, SliceNumber,:);
    redChannel=(reshape(redChannel,[DS.Dimension(1) DS.Dimension(3)])).*255;
    redChannel=imrotate(redChannel,-90);
    redChannel=(flipud(fliplr(redChannel)));
    greenChannel = DS.IPF001(2, :, SliceNumber,:);
    greenChannel=(reshape(greenChannel,[DS.Dimension(1) DS.Dimension(3)])).*255;
    greenChannel=imrotate(greenChannel,-90);
    greenChannel=(flipud(fliplr(greenChannel)));
    blueChannel = DS.IPF001(3, :, SliceNumber,:);
    blueChannel=(reshape(blueChannel,[DS.Dimension(1) DS.Dimension(3)])).*255;
    blueChannel=imrotate(blueChannel,-90);
    blueChannel=(flipud(fliplr(blueChannel)));
    Im_GB_rec_cl=joinchannels('RGB',redChannel, greenChannel, blueChannel);
    % add on Oct 1, 2020
    Im_GB_rec_bin=Im_GB_rec>0;
    Im_mask=~Im_GB_rec_bin;
    Im_mask=bwareaopen(Im_mask,round(size(Im_mask,1)*size(Im_mask,2)*0.01));
    Im_mask=Im_mask*255; % white background
    redChannel_mask=Im_mask;
    greenChannel_mask=Im_mask;
    blueChannel_mask=Im_mask;
    Im_mask_cl=joinchannels('RGB',redChannel_mask, greenChannel_mask, blueChannel_mask);
    Im_GB_rec_cl=Im_GB_rec_cl+Im_mask_cl;
elseif compare_YZ==1
    slice_rec=DS.GIDvol(SliceNumber,:,:);
    Im_GB_rec=(reshape(slice_rec,[DS.Dimension(2) DS.Dimension(3)]));
    Im_GB_rec=imrotate(Im_GB_rec,-90);
    Im_GB_rec=(flipud(fliplr(Im_GB_rec)));
    Dimension=[DS.Dimension(2) DS.Dimension(3)];
    redChannel = DS.IPF001(1, SliceNumber,:,:);
    redChannel=(reshape(redChannel,[DS.Dimension(2) DS.Dimension(3)])).*255;
    redChannel=imrotate(redChannel,-90);
    redChannel=(flipud(fliplr(redChannel)));
    greenChannel = DS.IPF001(2, SliceNumber,:,:);
    greenChannel=(reshape(greenChannel,[DS.Dimension(2) DS.Dimension(3)])).*255;
    greenChannel=imrotate(greenChannel,-90);
    greenChannel=(flipud(fliplr(greenChannel)));
    blueChannel = DS.IPF001(3,SliceNumber,:,:);
    blueChannel=(reshape(blueChannel,[DS.Dimension(2) DS.Dimension(3)])).*255;
    blueChannel=imrotate(blueChannel,-90);
    blueChannel=(flipud(fliplr(blueChannel)));
    Im_GB_rec_cl=joinchannels('RGB',redChannel, greenChannel, blueChannel);
    % add on Oct 1, 2020
    Im_GB_rec_bin=Im_GB_rec>0;
    Im_mask=~Im_GB_rec_bin;
    Im_mask=bwareaopen(Im_mask,round(size(Im_mask,1)*size(Im_mask,2)*0.01));
    Im_mask=Im_mask*255; % white background
    redChannel_mask=Im_mask;
    greenChannel_mask=Im_mask;
    blueChannel_mask=Im_mask;
    Im_mask_cl=joinchannels('RGB',redChannel_mask, greenChannel_mask, blueChannel_mask);
    Im_GB_rec_cl=Im_GB_rec_cl+Im_mask_cl;
end

