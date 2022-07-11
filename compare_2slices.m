% compare two slices from two different datasets
% DS_ref and DS_rec
% compare_XY, compare_XZ and compare_YZ define the slice in which plane
% grains: number of grains [1*2]
function [GB_overlay,DistPixel1]=compare_2slices(DS_ref,DS_rec, ...
    SliceNumber_ref,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ)

% % for testing
% compare_XY=1;
% compare_XZ=0;
% compare_YZ=0;
% SliceNumber_ref=40;
% SliceNumber_rec=40;

grains=[length(DS_ref.SeedID) length(DS_rec.SeedID)];
crop_edge_flag=0; % whether to remove the empty sample on the edge, 1: yes; 0: no

%%% change the default settings
setMTEXpref('xAxisDirection','east'); % default: 'north'
setMTEXpref('zAxisDirection','outOfPlane'); % same as default
setMTEXpref('bAxisDirection','north'); % default: 'east'
setMTEXpref('aAxisDirection',''); % same as default
setMTEXpref('FontSize',44); % default: 15
% define symmetries
cs = crystalSymmetry('cubic');
cK = ipfHSVKey(cs);
cK.inversePoleFigureDirection=zvector;

if compare_XY==1
    Slice{1}=DS_ref.GIDvol(:,:,SliceNumber_ref);
    Slice{2}=DS_rec.GIDvol(:,:,SliceNumber_rec);
    Im{1}=(reshape(Slice{1},[DS_ref.Dimension(1) DS_ref.Dimension(2)]));
    Im{2}=(reshape(Slice{2},[DS_rec.Dimension(1) DS_rec.Dimension(2)]));
    Dimension=[DS_ref.Dimension(1) DS_ref.Dimension(2)];
elseif compare_XZ==1
    Slice{1}=DS_ref.GIDvol(:,SliceNumber_ref,:);
    Slice{2}=DS_rec.GIDvol(:,SliceNumber_rec,:);
    Im{1}=(reshape(Slice{1},[DS_ref.Dimension(1) DS_ref.Dimension(3)]));
    Im{2}=(reshape(Slice{2},[DS_rec.Dimension(1) DS_rec.Dimension(3)]));
    Dimension=[DS_ref.Dimension(1) DS_ref.Dimension(3)];
elseif compare_YZ==1
    Slice{1}=DS_ref.GIDvol(SliceNumber_ref,:,:);
    Slice{2}=DS_rec.GIDvol(SliceNumber_rec,:,:);
    Im{1}=(reshape(Slice{1},[DS_ref.Dimension(2) DS_ref.Dimension(3)]));
    Im{2}=(reshape(Slice{2},[DS_rec.Dimension(2) DS_rec.Dimension(3)]));
    Dimension=[DS_ref.Dimension(2) DS_ref.Dimension(3)];
end
slice_ref=slice_show(DS_ref,SliceNumber_ref,compare_XY,compare_XZ,compare_YZ);
slice_rec=slice_show(DS_rec,SliceNumber_rec,compare_XY,compare_XZ,compare_YZ);

for j=1:length(Im)
    if crop_edge_flag==1
        Im_GB=Im{j};
        Im_GB=imrotate(Im_GB,-90);
        Im_GB=(flipud(fliplr(Im_GB)));
        if sum(all(Im_GB))~=Dimension(1) || sum(all(Im_GB))==0
            if sum(~all(Im_GB==0,1))==length(Im_GB(1,:))
                Index_dim(j,1)=1;
                Index_dim(j,2)=length(Im_GB(1,:));
            else
                for i=1:length(Im_GB(1,:))-1
                    if sum(Im_GB(:,i))==0 && sum(Im_GB(:,i+1))>0
                        Index_dim(j,1)=i+1;
                    end
                    if sum(Im_GB(:,i))>0 && sum(Im_GB(:,i+1))==0
                        Index_dim(j,2)=i;
                    end
                end
            end
            if Index_dim(j,1)==0
                Index_dim(j,1)=1;
            end
            Index_dim(j,:)=Index_dim(1,:); % This means the crop keeps the same for the reference image

            if sum(all(Im_GB))==0 && sum(~all(Im_GB==0,2))~=Dimension(2)
                NewDim=[Index_dim(j,2)-Index_dim(j,1)+1 sum(~all(Im_GB==0,2))];
                for k=1:length(Im_GB(:,1))-1
                    if all(Im_GB(k,:)==0) && ~all(Im_GB(k+1,:)==0)
                        Index_top(j,1)=k;
                    elseif ~all(Im_GB(k,:)==0) && all(Im_GB(k+1,:)==0)
                        Index_top(j,2)=k;
                    end
                end
                if Index_top(j,1)<1
                    Index_top(j,1)=1;
                    Index_top(j,2)=Index_top(j,2)+1;
                end
                Index_top(j,:)=Index_top(1,:);
                Im_GB=Im_GB(Index_top(j,1):Index_top(j,2),Index_dim(j,1):Index_dim(j,2));
%                     if all(Im_GB(end,:)==0)
%                         Im_GB=Im_GB(1:sum(~all(Im_GB==0,2)),Index_dim(j,1):Index_dim(j,2));
%                     else
%                         Im_GB=Im_GB(length(Im_GB(:,1))-sum(~all(Im_GB==0,2))+1:length(Im_GB(:,1)), ...
%                             Index_dim(j,1):Index_dim(j,2));
%                     end
            else
                NewDim=Dimension;
                Im_GB=Im_GB(:,Index_dim(j,1):Index_dim(j,2));
            end
        else
            NewDim=Dimension;
            Index_dim(j,1)=1;
            Index_dim(j,2)=length(Im_GB(1,:));
            Im_GB=Im_GB(:,Index_dim(1,1):Index_dim(1,2)); % add on Oct 17, 2019
        end
    else
        Im_GB=Im{j};
        Im_GB=imrotate(Im_GB,-90);
        Im_GB=(flipud(fliplr(Im_GB)));
        NewDim=Dimension;
        Index_dim(j,1)=1;
        Index_dim(j,2)=length(Im_GB(1,:));
    end
    Im_GB_rec{j}=Im_GB;
%         dipshow(Im_GB_rec{j});
end

% extract binary image for grain boundaries
for j=1:length(Im)
    Im_GB_rec_bin{j}=newim([length(Im_GB_rec{j}(1,:)),length(Im_GB_rec{j}(:,1))]);
    for i=1:grains(j)
        filter_index=[];
        FilterGrain=zeros(size(Im_GB_rec{j}));
        [filter_index(:,1) filter_index(:,2)]=ind2sub(size(Im_GB_rec{j}),(find(Im_GB_rec{j}==i)));
        if ~isempty(filter_index)
            for k=1:length(filter_index(:,1))
                FilterGrain(filter_index(k,1),filter_index(k,2))=1;
            end
            [Bgb,Lgb] = bwboundaries(FilterGrain,'noholes');
            for m=1:length(Bgb)
                for n=1:length(Bgb{m}(:,1))
                    Im_GB_rec_bin{j}(Bgb{m}(n,2)-1,Bgb{m}(n,1)-1)=1;
                end
            end
        end
        i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%% link color with orientation for the reference
%%%%%%%%%%%%%%%%%%%%%%% slice
GB_overlay=slice_ref;
clear Index_GB;
[Index_GB(:,1),Index_GB(:,2)]=(ind2sub(size(double(Im_GB_rec_bin{2})),(find(double(Im_GB_rec_bin{2})==1))));
for i=1:length(Index_GB(:,1))
    if sum(size(GB_overlay)==[3 1])==2
        GB_overlay(Index_GB(i,2)-1,Index_GB(i,1)-1)=[0 0 0]; % RGB color image
    else
        GB_overlay(Index_GB(i,2)-1,Index_GB(i,1)-1)=0; % grey-value image
    end
end
%     dipshow(Im_GB_rec{1}); % reference image
dipshow(GB_overlay); % overlay of the two

%%%%% quantify the difference between GB positions
% distance with respect to Im_ref
Dist=bwdist(double(Im_GB_rec_bin{1}),'euclidean'); % the GB pixels are nonzeros in the binary image
    
% method 1
filter_index=[];
[filter_index(:,1),filter_index(:,2)]=ind2sub([length(Im_GB_rec_bin{2}(1,:)) length(Im_GB_rec_bin{2}(:,1))], ...
    (find(Im_GB_rec_bin{2}==1)+1));
% find all the difference for the reconstructred structure
DistIm1=zeros(size(Dist));
for i=1:length(filter_index(:,1))
    DistIm1(filter_index(i,1),filter_index(i,2))=Dist(filter_index(i,1),filter_index(i,2));
end
DistPixel1=[];
for i=1:length(DistIm1(1,:))
    DistPixel1=[DistPixel1;DistIm1(:,i)];
end

% method 2
DistIm2=Dist.*double(Im_GB_rec_bin{2});
DistPixel2=[];
for i=1:length(DistIm2(1,:))
    DistPixel2=[DistPixel2;DistIm2(:,i)];
end

figure('Name','Statistics');
hh=histogram(DistPixel1,'BinWidth',1);
hh.Normalization = 'probability';
xlabel('Deviation (pixels)','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',16);
set(gca,'YScale','log');
box on;
title(strcat('Slice number',num2str(SliceNumber_ref)));
% export statistics of difference in GB position
% dlmwrite([strcat('compare2rec_',num2str(SliceNumber_ref)) '.txt'],[hh.BinEdges(1:end-1)' hh.Values'],'delimiter',' ');

