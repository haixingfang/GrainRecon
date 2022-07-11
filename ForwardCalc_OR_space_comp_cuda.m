% discretize the quaternion space
% perform forward simulation
% check the frequency of intersection
% July 12, 2021
% gpu cuda, April 8, 2022

function Output = ForwardCalc_OR_space_comp_cuda(pos,OR,S,B,Ahkl,nrhkl, ...
                        proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                        P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                        BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max)
% % for testing
% pos=pos_indexing;
hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';

nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
pos_cuda = single(pos);
q_cuda = single(reshape(OR.q',1,[]));
S_cuda = single(reshape(S',1,9));
B_cuda = single(reshape(B',1,9));
hkl_cuda = single(reshape(hkl,[1 size(hkl,1)*size(hkl,2)]));
proj_bin_bw_cuda = single(reshape(proj_bin_bw,1,[]));
rot_angles_cuda = single(rot_angles);
RotDet_cuda = single(reshape(RotDet',1,9));
param_cuda = single([nr_hkl nr_rot Lsam2sou Lsam2det P0y P0z dety00 detz00 ...
                    pixelysize pixelzsize detysize detzsize BeamStopY ...
                    BeamStopZ thetamax lambda_min lambda_max]);
if nr_rot * nr_hkl > 181*80
    error("error: the number of possible diffraction events exceeds the limit, \nplease increase the limit of the pre_allocated size for dis_simu_all in cuda_forward_comp.cu")
end

Output_cuda = cuda_forward_comp(pos_cuda, q_cuda, S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda);
tmp = gather(Output_cuda);
Output = zeros(size(OR.q,1),9);
for i=1:size(OR.q,1)
    Output(i,9) = i;
    Output(i,1:4) = tmp((i-1)*4+1:(i-1)*4+4);
    Output(i,5:8) = OR.q(i,:);
end
% Output = sortrows(Output,3,'descend');
gpuDevice(1); % reset gpu device

    


