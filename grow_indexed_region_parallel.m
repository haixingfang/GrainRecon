function [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region_parallel(I,Imask,DisMedian,Ninter,pos_indices,Spot_L,drop_off,maxD, ...
        indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max, ...
        Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z,RotAxisOffset, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian)
        
%
% I : input image - completeness map
% J : logical output image of region
% pos_indices: [x,y,z] the position of the seedpoint
% drop_off : completeness decrease amount, by default 0.02
% minD : minimum acceptable completeness weighted center change [pixel]
%
% Given a bounding box whose dimensions are calculated according to the spot size for checking the voxels to grow or not
% March 30, 2022

% % for testing
% I = Completeness;
% Imask = Mask;
% DisMedian = Dismedian;
% i=1;
% Spot_L = mean(Spots_pair{i}(:,21))+3*std(Spots_pair{i}(:,21)); % bounding area of the spot [pixel^2]
% Spot_L = sqrt(Spot_L/pi); % bounding half length of the spot [pixel]

% %%%
x = pos_indices(1);
y = pos_indices(2);
z = pos_indices(3);
I(x,y,z) = indexed_comp;
J = zeros(size(I)); % Output
J(x,y,z) = 2;
Completeness_out = I;
DisMedian_out = DisMedian;
Ninter_out = Ninter;

seed_comp = I(x,y,z); % The seeding completeness of the segmented region
completeness_sum = seed_comp;
completeness_sum_pos = seed_comp*[x,y,z];
completeness_center0 = completeness_sum_pos/completeness_sum; % weighted center

Box_L = ceil(Spot_L*(mean([pixelysize pixelzsize])) / ((Lsam2det/Lsam2sou)*VoxSize)); % half width of the bounding box for the grain [voxel]
% define the boundary of the box
for i = 1:3
    if pos_indices(i)-Box_L > 0 
        BoxDim(i,1) = pos_indices(i)-Box_L;
    else
        BoxDim(i,1) = 1;
    end
    if pos_indices(i)+Box_L < size(I,i) 
        BoxDim(i,2) = pos_indices(i)+Box_L;
    else
        BoxDim(i,2) = size(I,i);
    end
end
[x_grid,y_grid,z_grid] = meshgrid(BoxDim(1,1):BoxDim(1,2),BoxDim(2,1):BoxDim(2,2),BoxDim(3,1):BoxDim(3,2));
voxel_waiting_growth = [reshape(x_grid,1,[])' reshape(y_grid,1,[])' reshape(z_grid,1,[])'];
sprintf('Investigate the bounding volume [%d %d; %d %d; %d %d] for the seeding voxel [%d %d %d] and check %d voxels for growth.', ...
    reshape(BoxDim',1,[]),pos_indices,length(voxel_waiting_growth(:,1)))

neg_list = zeros(length(voxel_waiting_growth(:,1)),7);
if length(voxel_waiting_growth(:,1))>2000
    parfor i=1:length(voxel_waiting_growth(:,1))
        xn = voxel_waiting_growth(i,1);
        yn = voxel_waiting_growth(i,2);
        zn = voxel_waiting_growth(i,3);
        
        % Check if the voxel is in the mask
        if Imask(xn,yn,zn)>0
            pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
            if simap_data_flag==1
                pos(1)=-pos(1);
                pos(2)=-pos(2);
            end
            [Nr_simu,Nr_intersect,dis_median]=forward_comp(UU,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            if (Nr_intersect/Nr_simu > I(xn,yn,zn) && dis_median<=min([DisMedian(xn,yn,zn) maxDmedian])) ...
                    || (abs(Nr_intersect/Nr_simu-I(xn,yn,zn))<=0.01 && Nr_intersect>Ninter(xn,yn,zn))
%                 J(xn,yn,zn)=1;                
                if Nr_intersect/Nr_simu >= seed_comp*(1-drop_off) % accept this voxel
%                     J(xn,yn,zn)=2;
%                     Completeness_out(xn,yn,zn) = Nr_intersect/Nr_simu;
%                     DisMedian_out(xn,yn,zn) = dis_median;
%                     Ninter_out(xn,yn,zn) = Nr_intersect;
                    neg_list(i,:) = [xn yn zn Nr_intersect/Nr_simu dis_median Nr_intersect 1];
                else
                    neg_list(i,:) = [xn yn zn Nr_intersect/Nr_simu dis_median Nr_intersect 0];
                end
            end
        end
    end
else
    for i=1:length(voxel_waiting_growth(:,1))
        xn = voxel_waiting_growth(i,1);
        yn = voxel_waiting_growth(i,2);
        zn = voxel_waiting_growth(i,3);
        
        % Check if the voxel is in the mask
        if Imask(xn,yn,zn)>0
            pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
            if simap_data_flag==1
                pos(1)=-pos(1);
                pos(2)=-pos(2);
            end
            [Nr_simu,Nr_intersect,dis_median]=forward_comp(UU,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            if (Nr_intersect/Nr_simu > I(xn,yn,zn) && dis_median<=min([DisMedian(xn,yn,zn) maxDmedian])) ...
                    || (abs(Nr_intersect/Nr_simu-I(xn,yn,zn))<=0.01 && Nr_intersect>Ninter(xn,yn,zn))
%                 J(xn,yn,zn)=1;                
                if Nr_intersect/Nr_simu >= seed_comp*(1-drop_off) % accept this voxel
%                     J(xn,yn,zn)=2;
%                     Completeness_out(xn,yn,zn) = Nr_intersect/Nr_simu;
%                     DisMedian_out(xn,yn,zn) = dis_median;
%                     Ninter_out(xn,yn,zn) = Nr_intersect;
                    neg_list(i,:) = [xn yn zn Nr_intersect/Nr_simu dis_median Nr_intersect 1];
                else
                    neg_list(i,:) = [xn yn zn Nr_intersect/Nr_simu dis_median Nr_intersect 0];
                end
            end
        end
    end   
end

ind=find(neg_list(:,7)==1); % voxels to be accept
if ~isempty(ind)
    for j=1:length(ind)
        xn = neg_list(ind(j),1);
        yn = neg_list(ind(j),2);
        zn = neg_list(ind(j),3);
        J(xn,yn,zn)=2;
        Completeness_out(xn,yn,zn) = neg_list(ind(j),4);
        DisMedian_out(xn,yn,zn) = neg_list(ind(j),5);
        Ninter_out(xn,yn,zn) = neg_list(ind(j),6);
    end
end

% Return the grown volume as logical matrix
J=J>1;
[indexed_indices(:,1),indexed_indices(:,2),indexed_indices(:,3)]=ind2sub(size(J),find(J>0));

if length(indexed_indices(:,1))>1
    for j=1:length(indexed_indices(:,1))
        completeness_sum = completeness_sum + ...
            Completeness_out(indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3));
        completeness_sum_pos = completeness_sum_pos + ...
            Completeness_out(indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3))* ...
            [indexed_indices(j,1),indexed_indices(j,2),indexed_indices(j,3)];
    end
end

completeness_center = completeness_sum_pos/completeness_sum; % update weighted center
if sqrt(sum((completeness_center-completeness_center0).^2))>max([maxD maxD/10000*length(indexed_indices(:,1))]) % relative for pixel number > 10000
    center = round(completeness_center);
    replace_center = 1;
else
    center = completeness_center0;
    replace_center = 0;
end





