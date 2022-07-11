% thin down the binarized thresholded spots
% Feb 10, 2021
function im_spot_bin=Binary_segmentation_parallel(imdata,imbin,BackGround)

percent=0.1;
im1=double(imdata);
im_label=bwlabel(double(imbin),4);

im_spot=zeros(size(im1));
parfor j=1:max(max(im_label))
    ConnectedID=im_label==j;
    ConnectedComp=ConnectedID.*im1;
    thres_int=BackGround+percent*(max(max(ConnectedComp))-BackGround);
    AcceptComp=ConnectedComp>thres_int;
    im_spot=im_spot+AcceptComp;
    j;
end
im_spot_bin=im_spot>0;
im_spot_bin=bwareaopen(im_spot_bin,1); % remove object with only 1 pixel

