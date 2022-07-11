
function [newIm LocImax npeaks]=splitConnectedSpots(ConnectedComp,BackGround,percent)
% % for testing
% j=874;
% ConnectedID=im_label==j;
% ConnectedComp=ConnectedID.*im1;
% BackGround=0;

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
    stop_flag=0;
    factor=0.7;
    iter=1;
    while stop_flag~=1
        newIm=zeros(size(z));
        for i=1:npeaks
%             if (LocImax(i,3)-BackGround)<10
%                 BW=grayconnected(z,LocImax(i,1),LocImax(i,2),max([(LocImax(i,3)-BackGround)*factor 1]));
%             else
                BW=grayconnected(z,LocImax(i,1),LocImax(i,2),(LocImax(i,3)-BackGround)*factor);
%             end
            newIm=newIm+BW;
        end
        newIm=imclose(newIm,strel('square',5));
        newIm=newIm>0;
        newIm_label=bwlabel(newIm,8);
        s_msr=regionprops(newIm_label,'centroid');
%         length(s_msr)
        if (length(s_msr)>=npeaks-1 && npeaks>3) || (length(s_msr)>=2 && npeaks<=3) || factor<=0.4
            stop_flag=1;
        else
            factor=factor-0.025;
        end
        iter=iter+1;
    end
    
    newIm0=zeros(size(newIm));
    splitSpots_label=bwlabel(newIm,8);
    nspot=max(max(splitSpots_label));
    for k=1:max(max(splitSpots_label))
        SingleSpot=splitSpots_label==k;
        SingleSpot=imclose(SingleSpot,strel('disk',5)); % close holes
        if sum(sum(SingleSpot))>sum(sum(ConnectedComp>0))*0.025 % effective spot must be larger than a certain size
            newIm0=newIm0+SingleSpot;
        else
            nspot=nspot-1;
        end
    end
    newIm=newIm0>0;
end

plot_flag=0;
if plot_flag==1
    dipshow(z);
    dipshow(IMaxLoc);
    dipshow(newIm);
    [npeaks factor max(max(splitSpots_label)) nspot]
end
    


