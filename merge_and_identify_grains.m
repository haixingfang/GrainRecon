% merge regions and identify grains
function [DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_identify_grains(DS,mtex_avail,cs,min_misori,proj_bin_bw,Spots, ...
                    rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                    RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,minComp,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                    tomo_scale,VoxSize,RecVolumePixel,simap_data_flag)

[DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness(DS,mtex_avail, ...
    cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
iter_merge=1;
Ngrain_iter(iter_merge)=Ngrain;
iter_flag=1;  % iterative merging
if iter_flag==1
    stop_iter=0;
    while stop_iter~=1
        [DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness(DS_merge,mtex_avail, ...
            cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
        iter_merge=iter_merge+1;
        Ngrain_iter(iter_merge)=Ngrain;
        if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=0
            stop_iter=1;
        end
    end
end

if iter_flag==1
    stop_iter=0;
    while stop_iter~=1
         [DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness_mtex(DS_merge,1, ...
            cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
        iter_merge=iter_merge+1;
        Ngrain_iter(iter_merge)=Ngrain;
        if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=0
            stop_iter=1;
        end
    end
end
sprintf("Before revising: un-indexed fraction = %d / %d = %.3f",length(find(DS_merge.Mask==1))-length(find(DS_merge.GrainId>0)), ...
    length(find(DS_merge.Mask==1)),(length(find(DS_merge.Mask==1))- ...
    length(find(DS_merge.GrainId>0)))/length(find(DS_merge.Mask==1)))

sprintf("Use cpu for filling indexing ...")
[DS_merge]=revise_all_unindexed_voxel(DS_merge,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,minComp,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                tomo_scale,RecVolumePixel,simap_data_flag,0);
gpu_avail=0;
if gpuDeviceCount("available")
    if gpuDevice().TotalMemory/(1024*1024*1024)>6
        gpu_avail=1;
    end
end

if length(find(DS_merge.Mask==1))-length(find(DS_merge.GrainId>0))<=30000 || gpu_avail==1
    if length(find(DS_merge.Mask==1))-length(find(DS_merge.GrainId>0))<=2000
        sprintf("Use cpu for filling indexing ...")
        [DS_merge]=revise_all_unindexed_voxel(DS_merge,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,minComp,dety00,detz00,P0y,P0z, ...
                        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        tomo_scale,RecVolumePixel,simap_data_flag,1);
    else
        sprintf("Use gpu_cuda for filling indexing ...")
        [DS_merge] = revise_all_unindexed_voxel_cuda(DS_merge,proj_bin_bw,rot_angles,S,B,Ahkl,nrhkl, ...
                        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00, ...
                        P0y,P0z,pixelysize,pixelzsize,detysize,detzsize,BeamStopY,BeamStopZ, ...
                        RecVolumePixel,tomo_scale,VoxSize,simap_data_flag);
    end
end

stop_iter=0;
while stop_iter~=1
     [DS_merge,Ngrain,Nregion,Inherit_region_nr,CentroidComp]=merge_and_update_completeness_mtex(DS_merge,1, ...
        cs,min_misori,proj_bin_bw,Spots,rot_start,rot_step,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
        RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
        tomo_scale,VoxSize,RecVolumePixel,simap_data_flag);
    iter_merge=iter_merge+1;
    Ngrain_iter(iter_merge)=Ngrain;
    if Ngrain_iter(iter_merge-1)-Ngrain_iter(iter_merge)<=0
        stop_iter=1;
    end
end

sprintf("After revising: un-indexed fraction = %d / %d = %.3f",length(find(DS_merge.Mask==1))-length(find(DS_merge.GrainId>0)), ...
    length(find(DS_merge.Mask==1)),(length(find(DS_merge.Mask==1))- ...
    length(find(DS_merge.GrainId>0)))/length(find(DS_merge.Mask==1)))


