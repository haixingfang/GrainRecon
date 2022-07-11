% revise indexing for the unindexed single voxel by replacing with neighboring
% indexed information
% Jan 20, 2022
function [DS_merge]=revise_single_unindexed_voxel(DS_merge)

% % DS_merge_out=DS_merge;
% unindexed_indices=[];
% [unindexed_indices(:,1),unindexed_indices(:,2),unindexed_indices(:,3)]=ind2sub(size(DS_merge.GrainId), ...
%     find(DS_merge.GrainId==0 & DS_merge.Mask==1));

% this is faster
im_bin=DS_merge.GrainId==0 & DS_merge.Mask==1;
im_label=bwlabeln(im_bin,4);
im_msr=regionprops3(im_label,'all');
unindexed_indices=im_msr.Centroid(find(im_msr.Volume==1),:);
unindexed_indices=unindexed_indices(:,[2 1 3]);

neigb=[1 1 0;1 -1 0;1 1 1;1 0 1;1 -1 1;1 1 -1;1 0 -1;1 -1 -1; ...
    1 0 0;0 1 0;0 -1 0;0 1 1;0 0 1;0 -1 1;0 1 -1;0 0 -1;0 -1 -1; ...
    -1 1 0;-1 -1 0;-1 1 1;-1 0 1;-1 -1 1;-1 1 -1;-1 0 -1;-1 -1 -1;-1 0 0]; % 26 neighbors
Isizes = size(DS_merge.GrainId);
count_revise=0;
if ~isempty(unindexed_indices)
    for i=1:length(unindexed_indices(:,1))
        ind=unindexed_indices(i,:);
        Comp=[];
        GrainId=[];
        Rodrigues=[];
        EulerZXZ=[];
        IPF001=[];
        Dismedian=[];
        Icorr=[];
        VisitFlag=[];
        count0=0;
        count1=0;
        for j=1:size(neigb,1)
            % Calculate the neighbour coordinate
            xn = ind(1) + neigb(j,1);
            yn = ind(2) + neigb(j,2);
            zn = ind(3) + neigb(j,3);
            
            % Check if neighbour is inside or outside the image
            ins=(xn>=1)&&(yn>=1)&&(zn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2))&&(zn<=Isizes(3));
            if ins>0
                ind_neib=[xn yn zn];
                count0=count0+1;
                if (DS_merge.Completeness(xn,yn,zn)>0 && DS_merge.Mask(xn,yn,zn)==1) || ...
                        (DS_merge.Completeness(xn,yn,zn)==0 && DS_merge.Mask(xn,yn,zn)==0)
                    count1=count1+1;
                end
            end
        end
        if count1==count0 % it is an isolated unindexed voxel
%             DS_merge_out.Completeness(ind(1),ind(2),ind(3))=DS_merge.Completeness(ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.GrainId(ind(1),ind(2),ind(3))=DS_merge.GrainId(ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.Rodrigues(:,ind(1),ind(2),ind(3))=DS_merge.Rodrigues(:,ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.EulerZXZ(:,ind(1),ind(2),ind(3))=DS_merge.EulerZXZ(:,ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.IPF001(:,ind(1),ind(2),ind(3))=DS_merge.IPF001(:,ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.Dismedian(ind(1),ind(2),ind(3))=DS_merge.Dismedian(ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.Icorr(ind(1),ind(2),ind(3))=DS_merge.Icorr(ind_neib(1),ind_neib(2),ind_neib(3));
%             DS_merge_out.VisitFlag(ind(1),ind(2),ind(3))=DS_merge.VisitFlag(ind_neib(1),ind_neib(2),ind_neib(3));
            
            DS_merge.Completeness(ind(1),ind(2),ind(3))=DS_merge.Completeness(ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.GrainId(ind(1),ind(2),ind(3))=DS_merge.GrainId(ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.Rodrigues(:,ind(1),ind(2),ind(3))=DS_merge.Rodrigues(:,ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.EulerZXZ(:,ind(1),ind(2),ind(3))=DS_merge.EulerZXZ(:,ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.IPF001(:,ind(1),ind(2),ind(3))=DS_merge.IPF001(:,ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.Dismedian(ind(1),ind(2),ind(3))=DS_merge.Dismedian(ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.Icorr(ind(1),ind(2),ind(3))=DS_merge.Icorr(ind_neib(1),ind_neib(2),ind_neib(3));
            DS_merge.VisitFlag(ind(1),ind(2),ind(3))=DS_merge.VisitFlag(ind_neib(1),ind_neib(2),ind_neib(3));
            count_revise=count_revise+1;
        end
    end
end
sprintf('%d single isolated unindexed voxels have been corrected.',count_revise)
