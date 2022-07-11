% only for spots with more than 1 peak
% Feb 24, 2021

%DGG   Second derivative in the gradient-direction
%
% SYNOPSIS:
%  out = dgg(in,sigma)
%
% Computes:
%  g = gradient(in,sigma);
%  h = hessian(in,sigma);
%  out = (g'*h*g) / (g'*g);

function [newIm, LocImax, npeaks]=splitBydgg(ConnectedComp)
% %for testing
% j=150;
% ConnectedID=im_label==j;
% ConnectedComp=ConnectedID.*im1;

minDis=10; % minimum distance threshold for certifying the local maximum

z=double(ConnectedComp);
[XLocIMax,YLocIMax] = ind2sub(size(z),min(find(z == max(z,[],'all'))));
IMaxLoc=imregionalmax(z,8);
IMaxLoc_label=bwlabel(IMaxLoc,8);
s=regionprops(IMaxLoc_label,'centroid');
LocImax=[];
LocImax=[LocImax;XLocIMax YLocIMax z(XLocIMax,YLocIMax)];
AcceptedLocIMax=[XLocIMax YLocIMax];
IMaxLoc=zeros(size(z));
npeaks=1;
if length(s)>1 % more than one peak
    for i=1:length(s)
        dis=pdist([s(i).Centroid(2) s(i).Centroid(1);AcceptedLocIMax],'euclidean');
%         dis=sqrt((s(i).Centroid(1)-YLocIMax)^2+(s(i).Centroid(2)-XLocIMax)^2) % distance to the maximum
        if all(dis>=minDis)
            LocImax=[LocImax;floor(s(i).Centroid(2)) floor(s(i).Centroid(1)) z(floor(s(i).Centroid(2)),floor(s(i).Centroid(1)))];
            AcceptedLocIMax=[AcceptedLocIMax;s(i).Centroid(2) s(i).Centroid(1)];
            npeaks=npeaks+1;
        end
    end
end

maxPeaks=4; % detectable overlapped peaks maximum, 3 or 4
if npeaks>=maxPeaks 
    Z = linkage(AcceptedLocIMax,'centroid');
    c = cluster(Z,'Maxclust',maxPeaks);
%     scatter(AcceptedLocIMax(:,1),AcceptedLocIMax(:,2),10,c)
    for m=1:maxPeaks
        cindex(m)=max(find(c==m));
        IMaxLoc(floor(AcceptedLocIMax(cindex(m),1)),floor(AcceptedLocIMax(cindex(m),2)))=m;
    end
    LocImax=LocImax(cindex,:);
    AcceptedLocIMax=AcceptedLocIMax(cindex,:);
    
%     [LocImax,index]=sortrows(LocImax,3,'ascend');
%     AcceptedLocIMax=AcceptedLocIMax(index,:);
    npeaks=maxPeaks;
else
    IMaxLoc(floor(AcceptedLocIMax(1,1)),floor(AcceptedLocIMax(1,2)))=1;
end

if npeaks==1
%     thres_int=BackGround+percent*(max(max(z))-BackGround);
%     newIm=z>thres_int;
    newIm=zeros(size(z));
else
%     a=dip_image(z);
%     out=dgg(a,1);

%     % out = (g'*h*g) / (g'*g);
%     [gx,gy] = imgradientxy(z);
%     [gxx,gxy] = imgradientxy(gx);
%     [gyx,gyy] = imgradientxy(gy);   
    [gx, gy] = gradient(z); % this works best
    [gxx, gxy] = gradient(gx);
    [gyx, gyy] = gradient(gy);
    out=((gx.*gxx+gy.*gyx).*gx+(gx.*gxy+gy.*gyy).*gy)./(gx.^2+gy.^2);
    out(isnan(out))=0;
    
    out_bin=out<0;
    out_bin=double(out_bin);

    % removing small objects and then closing holes
    out_bin=bwareaopen(out_bin,5); % remove <5 pixels
    newIm=zeros(size(z));
    splitSpots_label=bwlabel(double(out_bin),8);
    nspot=max(max(splitSpots_label));
    for k=1:max(max(splitSpots_label))
        SingleSpot=splitSpots_label==k;
        SingleSpot=imclose(SingleSpot,strel('disk',5)); % close holes
        if sum(sum(SingleSpot))>sum(sum(ConnectedComp>0))*0.025 % effective spot must be larger than a certain size
            newIm=newIm+SingleSpot;
        else
            nspot=nspot-1;
        end
    end
    newIm=newIm>0;
end

plot_flag=0;
if plot_flag==1
    dipshow(ConnectedComp)
    dipshow(newIm);
    [npeaks max(max(splitSpots_label)) nspot]
end

