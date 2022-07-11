% convert pos_indexing [mm] to pos_indices
function pos_indices=pos2ind(i,pos_indexing,tomo_scale,VoxSize,RecVolumePixel,sz,simap_data_flag)
% sz=size(Completeness);

if simap_data_flag==1
    pos_indices(1)=(-pos_indexing(1)-tomo_scale.Center(1))/VoxSize+tomo_scale.Dimension(1)/2;
    pos_indices(2)=(-pos_indexing(2)-tomo_scale.Center(2))/VoxSize+tomo_scale.Dimension(2)/2;
    pos_indices(3)=(pos_indexing(3)-tomo_scale.Center(3))/VoxSize+tomo_scale.Dimension(3)/2;
else
    pos_indices=(pos_indexing-tomo_scale.Center')/VoxSize+tomo_scale.Dimension/2;
end
pos_indices=pos_indices-RecVolumePixel(:,1)'+1; % indices for the reconstructed box
% pos_indices=round(pos_indices);
pos_indices=ceil(pos_indices); % Dec 30, 2021

if pos_indices(1)>sz(1)
    pos_indices(1)=sz(1);
    warning('pos_indices %d may have exceeded the array bound, must not exceed %d',i,sz(1));
end
if pos_indices(2)>sz(2)
    pos_indices(2)=sz(2);
    warning('pos_indices %d may have exceeded the array bound, must not exceed %d',i,sz(2));
end
if pos_indices(3)>sz(3)
    pos_indices(3)=sz(3);
    warning('pos_indices %d may have exceeded the array bound, must not exceed %d',i,sz(3));
end

if pos_indices(1)<1
    pos_indices(1)=1;
    warning('pos_indices %d may have exceeded the array bound, must not smaller than 1',i);
end
if pos_indices(2)<1
    pos_indices(2)=1;
    warning('pos_indices %d may have exceeded the array bound, must not smaller than 1',i);
end
if pos_indices(3)<1
    pos_indices(3)=1;
    warning('pos_indices %d may have exceeded the array bound, must not smaller than 1',i);
end

