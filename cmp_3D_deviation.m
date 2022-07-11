% compute deviations for each pixel in 3D
% DS_ref and DS_rec
% compare_XY, compare_XZ and compare_YZ define the slice in which plane
% Output: DistIm - 3D deviation map
%         DistPixel - deviation for each voxel
%         dev_grain - deviation for each grain
%         pixel_dev_grain - deviation values for each pixel within one certain grain
function [DistIm,DistPixel,dev_grain,pixel_dev_grain]=cmp_3D_deviation(DS_ref,DS_rec,PairIndex, ...
                        SliceNumber,compare_XY,compare_XZ,compare_YZ)

% extract GB boundaries
Im_GB_ref=zeros(size(DS_ref.GIDvol));
for i=1:length(DS_ref.SeedID)
    id_vol=DS_ref.GIDvol==i;
    BW = edge3(id_vol,"approxcanny",0.2,[1 1 1]);
    ind=find(BW>0);
    Im_GB_ref(ind)=1;
%     i
end
Im_GB_rec=zeros(size(DS_rec.GIDvol));
for i=1:length(DS_rec.SeedID)
    id_vol=DS_rec.GIDvol==i;
    BW = edge3(id_vol,"approxcanny",0.2,[1 1 1]);
    ind=find(BW>0);
    Im_GB_rec(ind)=1;
%     i
end
% distance with respect to Im_ref
Dist=bwdist(double(Im_GB_ref),'euclidean'); % the GB pixels are nonzeros in the binary image
% assign deviation values for each pixels
DistIm=zeros(size(Dist));
DistIm(DS_ref.Mask==0)=NaN;

% find different pixels
dev_grain=zeros(length(PairIndex(:,1)),9);
for i=1:length(PairIndex(:,1))
    pixel_ind_ref=find(DS_ref.GIDvol==PairIndex(i,1));
    pixel_ind_rec=find(DS_rec.GIDvol==PairIndex(i,2));
    pixel_dev1=setdiff(pixel_ind_ref,pixel_ind_rec,'rows'); % (A, B) in A that is not in B
    pixel_dev2=setdiff(pixel_ind_rec,pixel_ind_ref,'rows');
   
    pixel_dev=[pixel_dev1;pixel_dev2];
    pixel_dev=unique(pixel_dev,'rows');
    if i==1 || i==length(PairIndex(:,1))
        sprintf('For grain %d, %d / %d = %.3f pixels are deviated from the reference data.',PairIndex(i,1),length(pixel_dev(:,1)),length(find(DS_ref.GIDvol==PairIndex(i,1))), ...
                            length(pixel_dev(:,1))/length(find(DS_ref.GIDvol==PairIndex(i,1))))
    end
    DistIm(pixel_dev)=Dist(pixel_dev);
    
    % for individual grains
    pixel_same=setdiff(pixel_ind_rec,pixel_dev2);
    pixel_dev_grain{i}=[zeros(length(pixel_same),1);Dist(pixel_dev2)];
    dev_grain(i,:)=[PairIndex(i,1:2) length(pixel_dev_grain{i}) length(find(pixel_dev_grain{i}>0)) ...
                length(find(pixel_dev_grain{i}>0))/length(pixel_dev_grain{i}) mean(pixel_dev_grain{i}) ...
                std(pixel_dev_grain{i}) mean(pixel_dev_grain{i}(pixel_dev_grain{i}>0)) ...
                std(pixel_dev_grain{i}(pixel_dev_grain{i}>0))];
end
DistPixel=DistIm((DistIm>=0));

% 2D visualization of middle slice
% SliceNumber=81;
% compare_XY=0;
% compare_XZ=1;
% compare_YZ=0;
if compare_XY==1
    im_slice=DistIm(:,:,SliceNumber);
    im_slice=reshape(im_slice,[DS_ref.Dimension(1) DS_ref.Dimension(2)]);
elseif compare_XZ==1
    im_slice=DistIm(:,SliceNumber,:);
    im_slice=reshape(im_slice,[DS_ref.Dimension(1) DS_ref.Dimension(3)]);
elseif compare_YZ==1
    im_slice=DistIm(SliceNumber,:,:);
    im_slice=reshape(im_slice,[DS_ref.Dimension(2) DS_ref.Dimension(3)]);
end
im_slice=imrotate(im_slice,-90);
im_slice=(flipud(fliplr(im_slice)));

figure;
imagesc(im_slice);
axis equal;
axis off;
box off;
colormap jet;
colorbar('FontSize',16);
set(gca,'visible','off');

