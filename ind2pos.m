% convert pos_indices to pos_indexing [mm]
function pos=ind2pos(pos_ind,tomo_scale,RecVolumePixel,simap_data_flag)
xn=pos_ind(1);
yn=pos_ind(2);
zn=pos_ind(3);
pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*tomo_scale.VoxSize(1)+tomo_scale.Center'; % [mm]
if simap_data_flag==1
    pos(1)=-pos(1);
    pos(2)=-pos(2);
end