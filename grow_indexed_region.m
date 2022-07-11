function [J,Completeness_out,DisMedian_out,Ninter_out,center,replace_center,indexed_indices]=grow_indexed_region(I,Imask,DisMedian,Ninter,pos_indices,drop_off,maxD, ...
        indexed_comp,UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max, ...
        Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian)
        
% This function performs "grow_indexed_region" in an image from a specified seedpoint (x,y,z)
% July 23, 2021
%
% I : input image - completeness map
% J : logical output image of region
% pos_indices: [x,y,z] the position of the seedpoint
% drop_off : completeness decrease amount, by default 0.02
% minD : minimum acceptable completeness weighted center change [pixel]
%
% The region iteratively grows by comparing all unallocated neighbouring pixels to the seeding completeness.
% The difference between a pixel's completeness value and the seeding value, 
% The pixel with the smallest change is allocated to the respective region. 
% This process stops when the completeness ratio between the seed and
% new pixel become larger than a certain treshold (drop_off)

% % for testing
% I = Completeness;
% Imask = Mask;

x = pos_indices(1);
y = pos_indices(2);
z = pos_indices(3);
I(x,y,z) = indexed_comp;
J = zeros(size(I)); % Output
Completeness_out = I;
DisMedian_out = DisMedian;
Ninter_out = Ninter;

Isizes = size(I); % Dimensions of input image
seed_comp = I(x,y,z); % The seeding completeness of the segmented region
reg_size = 0; % Number of pixels in region
completeness_sum = seed_comp;
completeness_sum_pos = seed_comp*[x,y,z];
completeness_center0 = completeness_sum_pos/completeness_sum; % weighted center
% Free memory to store neighbours of the (segmented) region
neg_free = 100;
neg_pos=0;
neg_list = zeros(neg_free,6);

pixdist=1; % initial difference in completeness value
% Neighbor locations (footprint)
% neigb=[-1 0 0; 1 0 0; 0 -1 0;0 1 0;0 0 1;0 0 -1]; % 6 neighbors
neigb=[1 1 0;1 -1 0;1 0 1;1 0 -1; ...
    1 0 0;0 1 0;0 -1 0;0 1 1;0 0 1;0 -1 1;0 1 -1;0 0 -1;0 -1 -1; ...
    -1 1 0;-1 -1 0;-1 0 1;-1 0 -1;-1 0 0]; % 18 neighbors
% neigb=[1 1 0;1 -1 0;1 1 1;1 0 1;1 -1 1;1 1 -1;1 0 -1;1 -1 -1; ...
%     1 0 0;0 1 0;0 -1 0;0 1 1;0 0 1;0 -1 1;0 1 -1;0 0 -1;0 -1 -1; ...
%     -1 1 0;-1 -1 0;-1 1 1;-1 0 1;-1 -1 1;-1 1 -1;-1 0 -1;-1 -1 -1;-1 0 0]; % 26 neighbors
stop_growing=0;
% Start regiogrowing until completeness value decrease by larger than a threshold value
% while (pixdist>1-drop_off && reg_size<sum(sum(sum(Imask>0))))
while stop_growing~=1
    % Add new neighbors pixels
    for j=1:size(neigb,1)
        % Calculate the neighbour coordinate
        xn = x + neigb(j,1);
        yn = y + neigb(j,2);
        zn = z + neigb(j,3);
        
        % Check if neighbour is inside or outside the image
        ins=(xn>=1)&&(yn>=1)&&(zn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2))&&(zn<=Isizes(3));
        
        % Add neighbor if inside and not already part of the segmented area
        if ins>0 && J(xn,yn,zn)==0 && Imask(xn,yn,zn)>0 %&& I(xn,yn,zn)<minComp
            pos=(([xn yn zn]+RecVolumePixel(:,1)'-1)-tomo_scale.Dimension/2).*VoxSize+tomo_scale.Center'; % [mm]
            if simap_data_flag==1
                pos(1)=-pos(1);
                pos(2)=-pos(2);
            end
            [Nr_simu,Nr_intersect,dis_median]=forward_comp(UU,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            if (Nr_intersect/Nr_simu > I(xn,yn,zn) && dis_median<=min([DisMedian(xn,yn,zn) maxDmedian])) ...
                    || (abs(Nr_intersect/Nr_simu-I(xn,yn,zn))<=0.01 && Nr_intersect>Ninter(xn,yn,zn))
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn zn Nr_intersect/Nr_simu dis_median Nr_intersect];
                J(xn,yn,zn)=1;
            end
        end
    end
    % Add a new block of free memory
    if(neg_pos+10>neg_free)
        neg_free=neg_free+100;
        neg_list(neg_pos+1:neg_free,:)=0;
    end
    
    if neg_pos==0 % no neighboring voxel to grow
        J(x,y,z)=2;
        stop_growing=1;
    else
        % Add pixel with smallest decrease of completeness
        dist = neg_list(1:neg_pos,4)/seed_comp;
        [pixdist, index] = max(dist);
        J(x,y,z)=2;
        reg_size=reg_size+1;

        if pixdist>1-drop_off && reg_size<sum(sum(sum(Imask>0)))-1
            % Calculate the new wighted center of the region
            completeness_sum = completeness_sum + neg_list(index,4);
            completeness_sum_pos = completeness_sum_pos + neg_list(index,4)*neg_list(index,1:3);
            Completeness_out(neg_list(index,1),neg_list(index,2),neg_list(index,3)) = neg_list(index,4); 
            DisMedian_out(neg_list(index,1),neg_list(index,2),neg_list(index,3)) = neg_list(index,5);
            Ninter_out(neg_list(index,1),neg_list(index,2),neg_list(index,3)) = neg_list(index,6);

            % Save the x and y coordinates of the pixel (for the neighbour add proccess)
            x = neg_list(index,1);
            y = neg_list(index,2);
            z = neg_list(index,3);

            % Remove the pixel from the neighbour (check) list
            neg_list(index,:)=neg_list(neg_pos,:);
            neg_pos=neg_pos-1;
        else
            stop_growing=1;
        end
    end
end

% Return the grown volume as logical matrix
J=J>1;
[indexed_indices(:,1),indexed_indices(:,2),indexed_indices(:,3)]=ind2sub(size(J),find(J>0));

completeness_center = completeness_sum_pos/completeness_sum; % update weighted center
if sqrt(sum((completeness_center-completeness_center0).^2))>max([maxD maxD/10000*length(indexed_indices(:,1))]) % relative for pixel number > 10000
    center = round(completeness_center);
    replace_center = 1;
else
    center = completeness_center0;
    replace_center = 0;
end





