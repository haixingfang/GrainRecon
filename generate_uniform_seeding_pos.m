% generate seeding positions from uniform grids on the mask
% replace the previous function generate_seeding_indexing_pos.m
% April 12, 2022
function pos_seed=generate_uniform_seeding_pos(iter,RecVolumePixel,tomo_scale,VoxSize,simap_data_flag)

grid_step_all = [50 45 40 35 30 25 20 15 10 6 3]; % 11 grid sizes
% grid_step_all = [50 40 35 30 25 20 17 13 9 6 3]; % 11 grid sizes,new

% nr=[];
% potential_voxels_all=[];
% for iter=1:10
grid_step = round(grid_step_all(iter));
dim = RecVolumePixel(:,2) - RecVolumePixel(:,1) +1;

if grid_step > min(dim)/3
    grid_step = max(round([min(dim)/3 grid_step/2]));
    if grid_step < 3
        grid_step = 3;
    end
end
for i = 1:3
    rand_number = 10*rand(1)-5;
    BoxDim(i,1) = RecVolumePixel(i,1)+round(grid_step/2+rand_number);
    if BoxDim(i,1) < RecVolumePixel(i,1)
        BoxDim(i,1) = RecVolumePixel(i,1);
    end
    BoxDim(i,2) = round(RecVolumePixel(i,2)-grid_step/2);
    if BoxDim(i,2) > RecVolumePixel(i,2)
        BoxDim(i,2) = RecVolumePixel(i,2);
    end
end
[x_grid,y_grid,z_grid] = meshgrid(BoxDim(1,1):grid_step:BoxDim(1,2),BoxDim(2,1):grid_step:BoxDim(2,2), ...
                                  BoxDim(3,1):grid_step:BoxDim(3,2));
potential_voxels = [reshape(x_grid,1,[])' reshape(y_grid,1,[])' reshape(z_grid,1,[])'];                              
% potential_voxels_all = [potential_voxels_all;potential_voxels];
% nr=[nr;length(potential_voxels(:,1))];
% end

seed_list = [];
for i=1:length(potential_voxels(:,1))
    if tomo_scale.Mask(potential_voxels(i,1),potential_voxels(i,2),potential_voxels(i,3))>0 ...
        && potential_voxels(i,1)>= RecVolumePixel(1,1) && potential_voxels(i,1)<=RecVolumePixel(1,2) ...
        && potential_voxels(i,2)>= RecVolumePixel(2,1) && potential_voxels(i,2)<=RecVolumePixel(2,2) ...
        && potential_voxels(i,3)>= RecVolumePixel(3,1) && potential_voxels(i,3)<=RecVolumePixel(3,2)
        seed_list=[seed_list;potential_voxels(i,:)];
    end
end
seed_list=seed_list(randperm(length(seed_list(:,1))),:); % make a random shuffle
seed_nr = length(seed_list(:,1));
sprintf('%d seeds are generated with an average spacing of %d pixels', seed_nr, grid_step)
if simap_data_flag==1
    pos_seed=(seed_list-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center'; % [mm] in lab system
    pos_seed(:,1)=-pos_seed(:,1); % applies to Simap data only
    pos_seed(:,2)=-pos_seed(:,2); % applies to Simap data only
else
    pos_seed=(seed_list-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center'; % [mm] in lab system
end
