% revise indexing for the single voxel by replacing with neighboring
% indexed information
% Dec 24, 2021
function [DS_merge]=revise_single_voxel(DS_new,DS_merge)

single_voxel_ID=find(DS_new.nVox==1);
neigb=[-1 0 0; 1 0 0; 0 -1 0;0 1 0;0 0 1;0 0 -1]; % 6 neighbors
Isizes = size(DS_new.GIDvol);
if ~isempty(single_voxel_ID)
    for i=single_voxel_ID'
        ind=[];
        [ind(1),ind(2),ind(3)]=ind2sub(size(DS_new.GIDvol),find(DS_new.GIDvol==i));
        Comp=[];
        GrainId=[];
        Rodrigues=[];
        EulerZXZ=[];
        IPF001=[];
        Dismedian=[];
        Icorr=[];
        VisitFlag=[];
        for j=1:size(neigb,1)
            % Calculate the neighbour coordinate
            xn = ind(1) + neigb(j,1);
            yn = ind(2) + neigb(j,2);
            zn = ind(3) + neigb(j,3);
            % Check if neighbour is inside or outside the image
            ins=(xn>=1)&&(yn>=1)&&(zn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2))&&(zn<=Isizes(3));
            if ins>0
                Comp=[Comp;DS_merge.Completeness(xn,yn,zn)];
                GrainId=[GrainId DS_merge.GrainId(xn,yn,zn)];
                Rodrigues=[Rodrigues;DS_merge.Rodrigues(:,xn,yn,zn)'];
                EulerZXZ=[EulerZXZ;DS_merge.EulerZXZ(:,xn,yn,zn)'];
                IPF001=[IPF001;DS_merge.IPF001(:,xn,yn,zn)'];
                Dismedian=[Dismedian;DS_merge.Dismedian(xn,yn,zn)];
                Icorr=[Icorr;DS_merge.Icorr(xn,yn,zn)];
                VisitFlag=[VisitFlag;DS_merge.VisitFlag(xn,yn,zn)];
            end
        end
        DS_merge.Completeness(ind(1),ind(2),ind(3))=mean(Comp);
        DS_merge.GrainId(ind(1),ind(2),ind(3))=GrainId(1);
        DS_merge.Rodrigues(:,ind(1),ind(2),ind(3))=Rodrigues(1,:);
        DS_merge.EulerZXZ(:,ind(1),ind(2),ind(3))=EulerZXZ(1,:);
        DS_merge.IPF001(:,ind(1),ind(2),ind(3))=IPF001(1,:);
        DS_merge.Dismedian(ind(1),ind(2),ind(3))=mean(Dismedian);
        DS_merge.Icorr(ind(1),ind(2),ind(3))=mean(Icorr);
        DS_merge.VisitFlag(ind(1),ind(2),ind(3))=VisitFlag(1);
    end
end

