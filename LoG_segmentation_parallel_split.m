
function im_spot_bin=LoG_segmentation_parallel_split(imdata,BackGround,sigma,percent,MinSpotSize)

im1=double(imdata);
% im1=medfilt2(im1);
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
im_label=bwlabel(im_LoG,8);

im_spot=zeros(size(im1));
algorithmKey=3; % [1 2 3]
parfor j=1:max(max(im_label))
    ConnectedID=im_label==j;
    ConnectedComp=ConnectedID.*im1;
    if sum(sum(ConnectedComp))>0
    if algorithmKey==1
        % algorithm 1, Feb 23, 2021
        [AcceptComp, LocImax, npeaks]=splitConnectedSpots(ConnectedComp,0,percent);
        if npeaks==1
            thres_int=percent*max(max(ConnectedComp));
            AcceptComp=ConnectedComp>thres_int;
        end
    elseif algorithmKey==2
        % algorithm 2, Feb 23, 2021
        [LocImax, npeaks]=FindNrOfPeaks(ConnectedComp);
        if npeaks==1
            thres_int=percent*max(max(ConnectedComp));
            AcceptComp=ConnectedComp>thres_int;
        else
            thres_int=percent*max(max(ConnectedComp));
            if thres_int>=0.95*min(LocImax(:,3))
                thres_int=0.95*min(LocImax(:,3));
            end
            AcceptComp=ConnectedComp>thres_int;
        end
    elseif algorithmKey==3
        % algorithm 3 based on dgg, Feb 23, 2021
        [AcceptComp, LocImax, npeaks]=splitBydgg(ConnectedComp);
        if npeaks==1
            thres_int=percent*max(max(ConnectedComp));
            AcceptComp=ConnectedComp>thres_int;
        end
    else
        thres_int=percent*max(max(ConnectedComp));
        AcceptComp=ConnectedComp>thres_int;        
    end
    im_spot=im_spot+AcceptComp;
    end
    
    j;
end
im_spot_bin=im_spot>0;
im_spot_bin=bwareaopen(im_spot_bin,MinSpotSize);
% split_label=label(im_spot_bin,1,MinSpotSize,size(imdata,1)*size(imdata,2))
% dipshow(imbin_bin-im_spot_bin)
