
% generate reconstructed slices from the tomography
% Cylinder shape: D*H
clear all;
dety0 = 80;	% beamcenter, y in pixel coordinatees
detz0 = 80;	% beamcenter, z in pixel coordinatees
detysize = 160;      % detector y size [pixels]
detzsize = 160;      % detector z size [pixels]
pixelysize=2.5e-3; % pixel size of the projection image [mm]
pixelzsize=2.5e-3; % pixel size of the projection image [mm]
direc='D:\Documents\Matlab\GrainRecon\Fe_100um_11_11_simu\tomo_slice';

D=400/1000; % [mm]
H=600/1000; % [mm]
nslice=round(H/pixelysize);
theta=0:0.005:2*pi;
SampleBoundary=[round(0.5*D/pixelysize.*cos(theta))+dety0;round(0.5*D/pixelzsize.*sin(theta))+detz0]';
SampleBoundary(SampleBoundary(:,1)==0,1)=1;
SampleBoundary(SampleBoundary(:,2)==0,2)=1;
for slicenumber=1:round(H/pixelysize)
    frame = zeros(detzsize,detysize);
    for i=min(SampleBoundary(:,1)):max(SampleBoundary(:,1))
        for j=min(SampleBoundary(:,2)):max(SampleBoundary(:,2))
            if sqrt((i-dety0)^2+(j-detz0)^2)<=round(0.5*D/pixelysize)
                frame(j,i)=(2^16-1)/2;
            end
        end
    end
    bgint=(2^16-1)/10;
    frame = frame + bgint*ones(detzsize,detysize); % Add background counts
    frame_image=uint16(frame);
    filename = sprintf('%s/tomo%0.4d.tif',direc,slicenumber-1); % Generate FILENAME of frame
    imwrite(frame_image,filename,'tif'); % Write out tiff file
    
%     filename = sprintf('%s/tomo%0.4d.raw',direc,slicenumber-1); % Generate FILENAME of frame
%     fid=fopen(filename,'w+');
%     fwrite(fid,frame_image,'uint16');
%     fclose(fid);
    
    slicenumber
end





