% find neiboring OR candidate
% April 20, 2022
function [id_neigb_all, id_neigb_ind] = find_neighbor_ORs(DS_in, wait_revising_voxels, search_radius)

nr_max = 30;
id_neigb_all = zeros(1,length(wait_revising_voxels(:,1))*nr_max); % store all the neighbors'ID
id_neigb_ind = zeros(length(wait_revising_voxels(:,1)),2);      % store the first and the last index of OR candidates for each voxel
count = 0;
for i=1:length(wait_revising_voxels(:,1))
    pos_indices = wait_revising_voxels(i,:);
    [ind,BoxDim] = find_ind(DS_in,i,pos_indices,search_radius);
%     ind=ind+BoxDim(:,1)'-1;
    if length(ind(:,1)) > nr_max
        sprintf('warning: the number of neighboring OR candidates %d is too many. Only use %d OR candidates', ...
            length(ind(:,1)),nr_max);
        nr_neigb = nr_max;
    else
        nr_neigb = length(ind(:,1));
    end

    id_neigb=zeros(1,nr_neigb);
    for m=1:nr_neigb
        id_neigb(m) = DS_in.GrainId(ind(m,1),ind(m,2),ind(m,3)); % neighbors ID
    end

    count = count + nr_neigb;
    id_neigb_all(count-nr_neigb+1:count) = id_neigb;
    id_neigb_ind(i,1) = count-nr_neigb+1;
    id_neigb_ind(i,2) = count;
end
id_neigb_all = id_neigb_all(1:id_neigb_ind(i,2)); % clean un-used space
% nr=id_neigb_ind(:,2)-id_neigb_ind(:,1);
% figure;hist(nr)
end


function [ind,BoxDim] = find_ind(DS_in,i,pos_indices,search_radius)

% define the boundary of the box
for j = 1:3
    if pos_indices(j)-search_radius > 0 
        BoxDim(j,1) = pos_indices(j)-search_radius;
    else
        BoxDim(j,1) = 1;
    end
    if pos_indices(j)+search_radius < size(DS_in.GrainId,j)
        BoxDim(j,2) = pos_indices(j)+search_radius;
    else
        BoxDim(j,2) = size(DS_in.GrainId,j);
    end
end
dim=BoxDim(:,2)-BoxDim(:,1)+1;
dim=dim';

clear ind;
[ind(:,1), ind(:,2), ind(:,3)] = ind2sub(dim,find(DS_in.GrainId(BoxDim(1,1):BoxDim(1,2), ...
    BoxDim(2,1):BoxDim(2,2),BoxDim(3,1):BoxDim(3,2))>0));
if isempty(ind)
    stop_flag=0;
    count=0;
    while stop_flag ~= 1
        count=count+1;
        % define the boundary of the box
        for j = 1:3
            if pos_indices(j)-(search_radius+count*2) > 0 
                BoxDim(j,1) = pos_indices(j)-(search_radius+count*2);
            else
                BoxDim(j,1) = 1;
            end
            if pos_indices(j)+(search_radius+count*2) < size(DS_in.GrainId,j)
                BoxDim(j,2) = pos_indices(j)+(search_radius+count*2);
            else
                BoxDim(j,2) = size(DS_in.GrainId,j);
            end
        end
        dim=BoxDim(:,2)-BoxDim(:,1)+1;
        dim=dim';
        clear ind;
        [ind(:,1), ind(:,2), ind(:,3)] = ind2sub(dim,find(DS_in.GrainId(BoxDim(1,1):BoxDim(1,2), ...
            BoxDim(2,1):BoxDim(2,2),BoxDim(3,1):BoxDim(3,2))>0));
        if ~isempty(ind)
            stop_flag=1;
        end
    end
    sprintf('Neighbor OR candidate is only found for voxel no. %d after increasing search_radius from %.1f to %.1f pixels', ...
        i,search_radius,search_radius+count*2)
end

% make the ind for neighboring grain IDs unique
ind=ind+BoxDim(:,1)'-1;
grainID_all=DS_in.GrainId(sub2ind(size(DS_in.GrainId),ind(:,1),ind(:,2),ind(:,3)));
[~,ia]=unique(grainID_all);
ind_eff=ind(ia,:);
ind=ind_eff;

end


