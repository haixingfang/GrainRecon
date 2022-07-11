% Calculate the completeness for all voxels which has been indexed with an orientation
% Suitable to check the completeness map for grain structure obtained from other techniques, e.g. synchrotron DCT
% July 8, 2022
function [DS_out] = Forward_calc_comp_allvoxel_fun(DS_in,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00, ...
                        P0y,P0z,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxD,maxDmedian)

DS_out = DS_in;
% initialize
Mask=DS_in.Mask;
Completeness=zeros(size(DS_in.GIDvol));
Dismedian=zeros(size(DS_in.GIDvol))+100;
Ninter=zeros(size(DS_in.GIDvol));

% indices for un-indexed voxels
for i = 1:length(DS_in.SeedID)
    un_indexed_voxels=[];
    if DS_in.SeedID(i) > 0
        [un_indexed_voxels(:,1), un_indexed_voxels(:,2), un_indexed_voxels(:,3)]=ind2sub(size(DS_in.GIDvol), ...
                find(DS_in.GIDvol==DS_in.SeedID(i)));
        if ~isempty(un_indexed_voxels)
            fprintf('Grain %d with %d voxels for computing completeness...\n',DS_in.SeedID(i),length(un_indexed_voxels(:,1)))
            UU=euler2u(DS_in.EulerZXZ(i,1)*pi/180,DS_in.EulerZXZ(i,2)*pi/180,DS_in.EulerZXZ(i,3)*pi/180);
            % call gpu function
            [Completeness_out,DisMedian_out,Ninter_out,Comp_max]=Forward_calc_comp_allvoxel_cuda(Completeness,Mask,Dismedian,Ninter,un_indexed_voxels,maxD, ...
                                UU,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                                RotDet,thetamax,lambda_min,lambda_max, ...
                                Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                                pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                                RecVolumePixel,tomo_scale,VoxSize,simap_data_flag,maxDmedian);
            % update
            Completeness=Completeness_out;
            Dismedian=DisMedian_out;
            Ninter=Ninter_out;
            DS_out.SeedComp(i)=Comp_max;
        else
            fprintf('find 0 voxels for grain no. %d\n',DS_in.SeedID(i))
        end
    end
end
fprintf('Done with %d voxels for %d grains.\n',length(find(DS_in.GIDvol>0)),length(find(DS_in.SeedID>0)))
DS_out.CompVol=Completeness;
DS_out.Dismedian=Dismedian;
DS_out.Ninter=Ninter_out;

