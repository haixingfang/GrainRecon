% remove spots from projections
% return new proj, proj_bin, proj_bin_bw
function [proj_out,proj_bin_out,proj_bin_bw_out]=remove_spots_from_proj(proj,proj_bin,SpotsForIndex)

for i=1:length(proj_bin)
    im=bwareaopen(proj_bin{i},3);
    im_label=bwlabel(im,8);
    im_measure=regionprops(im_label,'Area','Centroid','BoundingBox','EquivDiameter');
    spot_centroid=zeros(length(im_measure),2);
    for j=1:length(im_measure)
        spot_centroid(j,:)=im_measure(j).Centroid;
    end
    proj_new=zeros(size(proj{i}));
    spot_count=0;
    for j=1:length(SpotsForIndex{i})
       dis=sqrt(sum((repmat(SpotsForIndex{i}(j,9:10),length(spot_centroid(:,1)),1)-spot_centroid).^2,2));
       [dis_min,ind_spot]=min(dis);
       if dis_min<3
           spot_select=im_label==ind_spot;
           proj_new=proj_new+spot_select;
           spot_count=spot_count+1;
       end
    end
    proj_out{i}=proj{i}.*uint16(proj_new);
    proj_bin_out{i}=double(proj_out{i}>0);
    fprintf('Retrieved %d / %d spots for projection %d ...\n',spot_count,length(im_measure),i);
end
for i=1:length(proj_bin)
    [proj_bin_bw_out(:,:,i),~] = bwdist(double(proj_bin_out{i}));
end
