%%% post data analys
pos_seed_all_voxel=[];
for i=1:length(pos_seed_all(:,1))
    pos_indexing=pos_seed_all(i,:);
    if simap_data_flag==1
        pos_indices(1)=(-pos_indexing(1)-tomo_scale.Center(1))/VoxSize+tomo_scale.Dimension(1)/2;
        pos_indices(2)=(-pos_indexing(2)-tomo_scale.Center(2))/VoxSize+tomo_scale.Dimension(2)/2;
        pos_indices(3)=(pos_indexing(3)-tomo_scale.Center(3))/VoxSize+tomo_scale.Dimension(3)/2;
    else
        pos_indices=(pos_indexing-tomo_scale.Center')/VoxSize+tomo_scale.Dimension/2;
    end
    pos_indices=pos_indices-RecVolumePixel(:,1)'+1; % indices for the reconstructed box
    pos_indices=round(pos_indices);
    pos_seed_all_voxel=[pos_seed_all_voxel;pos_indices];
end

figure;
hold all;
plot3(pos_seed_all_voxel(1:seed_nr,1),pos_seed_all_voxel(1:seed_nr,2), ...
    pos_seed_all_voxel(1:seed_nr,3),'r.','MarkerSize',12,'LineWidth',1.5);
plot3(pos_seed_all_voxel(1*seed_nr+1:2*seed_nr,1),pos_seed_all_voxel(1*seed_nr+1:2*seed_nr,2), ...
    pos_seed_all_voxel(1*seed_nr+1:2*seed_nr,3),'bo','MarkerSize',6,'LineWidth',1.5);
plot3(pos_seed_all_voxel(2*seed_nr+1:end,1),pos_seed_all_voxel(2*seed_nr+1:end,2), ...
    pos_seed_all_voxel(2*seed_nr+1:end,3),'kx','MarkerSize',8,'LineWidth',1.5);
Mask = permute(DS.Mask,[2 1 3]);
ind=find(Mask==1);
AlphaValue=zeros(size(Mask));
AlphaValue(ind)=repmat(0.05,1,length(ind));
h1 = vol3d('cdata',Mask,'texture','3D','Alpha',AlphaValue);
view(3);
xlabel('x (pixel)');
ylabel('y (pixel)');
zlabel('z (pixel)');
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-27);
set(get(gca,'zlabel'),'rotation',90);
box on;
legend('iter 1','iter 2','iter 3');
set(gca,'FontSize',16);
set(gca,'LineWidth',1.5);
print(fullfile(FileFolder,'seeding_iter3'),'-dtiff','-r300');

i=3;
figure;
% h1 = vol3d('cdata',DS_record{i}.Completeness,'texture','3D');
h1 = vol3d('cdata',permute(DS_record{i}.Completeness,[2 1 3]),'texture','3D');
xlim([0 dim(1)]);
ylim([0 dim(2)]);
zlim([0 dim(3)]);
view(3);
% alphamap('vup');
colormap jet;
xlabel('x (pixel)');
ylabel('y (pixel)');
zlabel('z (pixel)');
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-27);
set(get(gca,'zlabel'),'rotation',90);
colorbar('FontSize',16);
set(gca,'FontSize',16);
set(gca,'LineWidth',1.5);
set(gca,'visible','on');
print(fullfile(FileFolder,['map_seed200_iter' num2str(i)]),'-dtiff','-r300');

figure;
hold all;
plot(hittedSpots_pair(:,13),hittedSpots_pair(:,14),'ro');
plot(hittedSpots_pair(:,19),hittedSpots_pair(:,20),'bx');
legend('Forward simu','Experimental');
xlim([0 detysize]);
ylim([0 detzsize]);
set(gca,'FontSize',16);
set(gca,'LineWidth',1.5);
box on;
print(fullfile(FileFolder,'centroids_fit'),'-dtiff','-r300');




