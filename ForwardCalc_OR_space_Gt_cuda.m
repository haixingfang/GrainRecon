% discretize the quaternion space
% perform forward simulation
% check the matched number of G vectors
% Sep 24, 2021
% gpu cuda, April 8, 2022

function Output = ForwardCalc_OR_space_Gt_cuda(pos,OR,S,B,Ahkl,nrhkl,hkl_family_square,d_possible, ...
                        Spots,rot_angles,RotDet,rot_start,rot_step,Lsam2sou,Lsam2det, ...
                        P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize)
% % for testing
% pos=pos_indexing;
hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';

hklnumber = length(d_possible);
nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
pos_cuda = single(pos);
q_cuda = single(reshape(OR.q',1,[]));
S_cuda = single(reshape(S',1,9));
B_cuda = single(reshape(B',1,9));
hkl_cuda = single(reshape(hkl,[1 size(hkl,1)*size(hkl,2)]));
Spots_cuda = [];
nr_spots = 0;
for i=1:length(Spots)
    Spots_cuda = [Spots_cuda;Spots{i}];
    nr_spots = nr_spots + length(Spots{i}(:,1));
end
Spots_cuda = single(reshape(Spots_cuda',1,[]));

rot_angles_cuda = single(rot_angles);
RotDet_cuda = single(reshape(RotDet',1,9));
hkl_family_square_cuda = single(hkl_family_square);
d_possible_cuda = single(d_possible);
param_cuda = single([hklnumber nr_hkl nr_rot nr_spots rot_start rot_step Lsam2sou Lsam2det P0y P0z ...
                    dety00 detz00 pixelysize pixelzsize detysize detzsize]);

% mexcuda cuda_forward_Gt.cu
% clear Output_cuda tmp Output;
Output_cuda = cuda_forward_Gt(pos_cuda, q_cuda, S_cuda, B_cuda, hkl_cuda, Spots_cuda, rot_angles_cuda, ...
                    RotDet_cuda, hkl_family_square_cuda, d_possible_cuda, param_cuda);
tmp = gather(Output_cuda);
Output = zeros(size(OR.q,1),8);
for i=1:size(OR.q,1)
    Output(i,1) = tmp((i-1)*4+1)+1;
    Output(i,2:4) = tmp((i-1)*4+2:(i-1)*4+4);
    Output(i,5:8) = OR.q(i,:);
end
% Output = sortrows(Output,2,'descend');
gpuDevice(1); % reset gpu device
    


