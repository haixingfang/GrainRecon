% find the indices for the voxels which have low completeness and small
% grain size
% Oct 13, 2022
function doubt_indexed_voxels = find_doubt_indexed_voxels(DS_in,minComp,small_grain_only)
doubt_indexed_voxels=[];
nVox_min=5;
minComp_thres=min([minComp*1.15 0.52]);

if small_grain_only~=1
    doubt_indexed_GrainId=DS_in.SeedID(DS_in.nVox<=nVox_min | DS_in.SeedComp<=minComp_thres);
else
    doubt_indexed_GrainId=DS_in.SeedID(DS_in.nVox<=nVox_min);
end

if ~isempty(doubt_indexed_GrainId)
    for i=doubt_indexed_GrainId
        ind=[];
        [ind(:,1),ind(:,2),ind(:,3)]=ind2sub(size(DS_in.GrainId),find(DS_in.GrainId==i));
        doubt_indexed_voxels=[doubt_indexed_voxels;ind];
    end
    fprintf('find %d doubtful grains (<= %d voxels or seedcomp <= %.3f): %d voxels\n', ...
        length(doubt_indexed_GrainId),nVox_min,minComp_thres,length(doubt_indexed_voxels(:,1)))
else
    fprintf('find 0 doubtful grains (<= %d voxels)\n',nVox_min)
end