% generate seeding positions for indexing
% new version
% April 23, 2022
function pos_seed=generate_seeding_indexing_pos(iter,allVoxel_indices,tomo_scale,VoxSize,simap_data_flag)

% seed_vol_min_all = [60 50 40 35 30 25 20 15 10 6 3].^3; % 11 sizes
% seed_vol_min_all = [45 40 35 30 25 22 19 16 13 10 6].^3; % 11 sizes,new
% seed_vol_min_all = [30 25 20 18 16 14 12 10 8 6 3].^3; % 11 sizes,new

seed_vol_min_all = [45 40 35 30 25 22 19 16 13 10 6 3].^3; % 12 sizes,new
if iter <= length(seed_vol_min_all)
    seed_vol_min = seed_vol_min_all(iter);
else
    seed_vol_min = seed_vol_min_all(end);
end
total_voxel=length(allVoxel_indices(:,1));
seed_nr=round(total_voxel/seed_vol_min);

if seed_nr<1
    seed_nr=round(total_voxel/(2*2*2));
elseif seed_nr>10000
    seed_nr=10000;
end
sprintf('%d seeds are suggested. Generating the seeds now ...',seed_nr)
dseed=(3*seed_vol_min/(4*pi))^(1/3)*2; % [pixel]
dmin=max([0.25*dseed 3]);          % [pixel]

allVoxel_indices=allVoxel_indices(randperm(length(allVoxel_indices(:,1))),:); % make a random shuffle
% sampling part of all voxels to speed up
if seed_nr<1000
    sampling_nr = min([round(seed_nr) length(allVoxel_indices(:,1))]);
elseif seed_nr<5000
    sampling_nr = min([round(seed_nr/3) length(allVoxel_indices(:,1))]);
elseif seed_nr<8000
    sampling_nr = min([round(seed_nr/6) length(allVoxel_indices(:,1))]);
else
    sampling_nr = min([round(seed_nr/50) length(allVoxel_indices(:,1))]);
end
allVoxel_indices_sample=allVoxel_indices(1:sampling_nr,:);

seed_list=allVoxel_indices_sample(1,:);
allVoxel_indices_pool=setdiff(allVoxel_indices_sample,seed_list,'rows');
allVoxel_indices_pool=allVoxel_indices_pool(randperm(length(allVoxel_indices_pool(:,1))),:); % make a random shuffle
count1=0;
count2=0;
seed_list_temp=[];
cutoff_nr=4000;
i=1;
while i<=seed_nr-1
   seed_list0=allVoxel_indices_pool(1,:);
   comb=[seed_list;seed_list0];
   distance=pdist(comb,'euclidean');
   if min(distance)>dmin
      seed_list=[seed_list;seed_list0];
      if i<cutoff_nr
          allVoxel_indices_pool=setdiff(allVoxel_indices_sample,seed_list,'rows');
      else
          allVoxel_indices_pool=setdiff(allVoxel_indices_sample,[seed_list_temp;seed_list],'rows');
      end
      if length(allVoxel_indices_pool(:,1))<50
          count1=count1+1;
          if sampling_nr*(count1+1)<length(allVoxel_indices(:,1)) && (1+count1*sampling_nr)<length(allVoxel_indices(:,1))
	          allVoxel_indices_sample=allVoxel_indices(1+count1*sampling_nr:sampling_nr*(count1+1),:);
          else
              sprintf('But no more available potential positions to be seeded.')
	          break;
          end
          if i<cutoff_nr
              allVoxel_indices_pool=setdiff(allVoxel_indices_sample,seed_list,'rows');
          else
              allVoxel_indices_pool=setdiff(allVoxel_indices_sample,[seed_list_temp;seed_list],'rows');
          end
      end
      i=i+1;
      if mod(i,cutoff_nr)==0 || i==seed_nr-1
        sprintf('Seed %d has been generated ...',i+1)
      end
   else
       count2=count2+1;
   end
   if count2>30
       count1=count1+1;
       if sampling_nr*(count1+1)<length(allVoxel_indices(:,1)) && (1+count1*sampling_nr)<length(allVoxel_indices(:,1))
           allVoxel_indices_sample=allVoxel_indices(1+count1*sampling_nr:sampling_nr*(count1+1),:);
       else
           sprintf('But no more available potential positions to be seeded.')
           break;
       end
       if i<cutoff_nr
           allVoxel_indices_pool=setdiff(allVoxel_indices_sample,seed_list,'rows');
       else
           allVoxel_indices_pool=setdiff(allVoxel_indices_sample,[seed_list_temp;seed_list],'rows');
       end
       count2=0;
   end
   allVoxel_indices_pool=allVoxel_indices_pool(randperm(length(allVoxel_indices_pool(:,1))),:);
   if mod(i,cutoff_nr)==0
       seed_list_temp=[seed_list_temp;seed_list(1:end-1,:)];
       seed_list=seed_list(end,:);
   end
end
if i>=cutoff_nr
    seed_list=[seed_list_temp;seed_list];
%     seed_list0=unique(seed_list,"rows");
end
% min(pdist(seed_list,'euclidean'))

sprintf('%d seeds are generated with a cubic dimension of %.1f pixels and dmin = %.1f pixels', ...
    length(seed_list(:,1)), seed_vol_min^(1/3),dmin)
if simap_data_flag==1
    pos_seed=(seed_list-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center'; % [mm] in lab system
    pos_seed(:,1)=-pos_seed(:,1); % applies to Simap data only
    pos_seed(:,2)=-pos_seed(:,2); % applies to Simap data only
else
    pos_seed=(seed_list-tomo_scale.Dimension/2)*VoxSize+tomo_scale.Center'; % [mm] in lab system
end

