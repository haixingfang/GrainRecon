
function [LocImax npeaks]=FindNrOfPeaks(ConnectedComp)

minDis=10; % minimum distance threshold for certifying the local maximum

z=double(ConnectedComp);
[XLocIMax,YLocIMax] = ind2sub(size(z),min(find(z == max(z,[],'all'))));
IMaxLoc=imregionalmax(z,8);
IMaxLoc_label=bwlabel(IMaxLoc,8);
s=regionprops(IMaxLoc_label,'centroid');
LocImax=[];
LocImax=[LocImax;XLocIMax YLocIMax z(XLocIMax,YLocIMax)];
AcceptedLocIMax=[XLocIMax YLocIMax];
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

if npeaks>=3 % detect 3 overlapped peaks maximum
    Z = linkage(AcceptedLocIMax,'centroid');
    c = cluster(Z,'Maxclust',3);
%     scatter(AcceptedLocIMax(:,1),AcceptedLocIMax(:,2),10,c)
    cindex(1)=max(find(c==1));
    cindex(2)=max(find(c==2));
    cindex(3)=max(find(c==3));
    LocImax=LocImax(cindex,:);
    AcceptedLocIMax=AcceptedLocIMax(cindex,:);
    
%     [LocImax,index]=sortrows(LocImax,3,'ascend');
%     AcceptedLocIMax=AcceptedLocIMax(index,:);
    npeaks=3;
end


