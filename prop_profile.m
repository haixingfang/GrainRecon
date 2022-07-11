% calculate the average and std along the distance to the center

function profile=prop_profile(prop,DS,grainno,id)

% prop=permute(prop,[2 1 3]);
center_vox=round(mean(id));
sz=size(prop);
R_min=0;
R_max=(round(DS.nVox(grainno).^(1/3))+1)/2;
% R_max=12;
R_thick=1; % [pixel]
R_center=R_min+R_thick/2:R_thick:R_max-R_thick/2;
[x, y, z] = meshgrid(1:sz(1), 1:sz(2), 1:sz(3));
distanceImage = sqrt((x-center_vox(1)).^2 + (y-center_vox(2)).^2 + (z-center_vox(3)).^2);
distanceImage=permute(distanceImage,[2 1 3]);
for j=1:length(R_center)
    clear indices;
   [indices(:,1),indices(:,2),indices(:,3)]=ind2sub(size(distanceImage), ...
       find(distanceImage>=R_center(j)-R_thick/2 & distanceImage<=R_center(j)+R_thick/2));
    im=zeros(sz(1),sz(2),sz(3));
    for k=1:length(indices(:,1))
        im(indices(k,1),indices(k,2),indices(k,3))=1;
    end
    im_select=double(im.*prop);
    im_nonzero=im_select(im_select>0);
    profile(j,1)=R_center(j);
    profile(j,2)=mean(im_nonzero);
    profile(j,3)=std(im_nonzero);
end


