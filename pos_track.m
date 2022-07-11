
clear all;
% close all;

FileFolder='D:\Documents\Matlab\GrainRecon\virtual_Fe_100um_6grains_rec\test_index_compete';
for i=1:10
    filename=['pos' num2str(i) '.mat'];
    load(fullfile(FileFolder,filename));
    pos{i}=pos_seed;
    pos_new{i}=pos_seed_new;
    
%     figure;
%     hold all;
%     % plot3(pos_all(:,1),pos_all(:,2),pos_all(:,3),'r.');
%     plot3(pos_seed(:,1),pos_seed(:,2),pos_seed(:,3),'ro');
%     plot3(pos_seed_new(:,1),pos_seed_new(:,2),pos_seed_new(:,3),'bx');
%     view(3);
end

figure;
hold all;
for i=1:length(pos)
    if ~isempty(pos{i})
        plot3(pos{i}(:,1),pos{i}(:,2),pos{i}(:,3),'o');
    end
end
view(3);

figure;
hold all;
for i=1:length(pos_new)
    if ~isempty(pos_new{i})
        plot3(pos_new{i}(:,1),pos_new{i}(:,2),pos_new{i}(:,3),'o');
    end
end
view(3);

