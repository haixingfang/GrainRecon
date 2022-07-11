
function im_spot_bin=LoG_segmentation_parallel(imdata,BackGround,sigma,percent,MinSpotSize)

im1=double(imdata);
im1=im1-BackGround; % step 1: subtract a defined background value
im1(im1<0)=0;
kernelSize=ceil(8*sqrt(2)*sigma);
G=fspecial('gaussian',[1,2*kernelSize+1],sigma);
d2G=G .*((-kernelSize:kernelSize).^2-sigma^2)/(sigma^4);
dxx = conv2(d2G, G, im1, 'same');
dyy = conv2(G, d2G, im1, 'same');
im2 = dxx + dyy;
im2 = -im2; 
im_LoG=im2>0; % step 2: apply LoG operator
im_LoG=bwareaopen(im_LoG,MinSpotSize);
% im_label=bwlabel(im_LoG,8);
im_label=bwlabel(im_LoG,4);

im_spot=zeros(size(im1));
parfor j=1:max(max(im_label))
    ConnectedID=im_label==j;
    ConnectedComp=ConnectedID.*im1;
    thres_int=percent*max(max(ConnectedComp));
    AcceptComp=ConnectedComp>thres_int;
    im_spot=im_spot+AcceptComp;
    j;
end
im_spot_bin=im_spot>0;
im_spot_bin=bwareaopen(im_spot_bin,MinSpotSize);

