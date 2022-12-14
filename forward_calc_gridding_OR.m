% forward calculate completeness for gridding orientation around the each of input orientation
% Oct 28, 2022

function [Output,Output_max] = forward_calc_gridding_OR(index_seeds,nr,gridsize,pos_cuda,S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda)

for j=1:length(index_seeds(:,1))
    OR_local=ori_local_sampling(index_seeds(j,5:8),nr,gridsize);
    if j==1
        q_all=zeros(OR_local.len*length(index_seeds(:,1)),5);
    end
    q_all((j-1)*length(OR_local.q(:,1))+1:j*length(OR_local.q(:,1)),1:4)=OR_local.q;
%     q_all((j-1)*length(OR_local.q(:,1))+1:j*length(OR_local.q(:,1)),5)=j;
    q_all((j-1)*length(OR_local.q(:,1))+1:j*length(OR_local.q(:,1)),5)=index_seeds(j,9);
end
q_cuda = single(reshape(q_all(:,1:4)',1,[]));
Output_cuda = cuda_forward_comp(pos_cuda, q_cuda, S_cuda, B_cuda, hkl_cuda, proj_bin_bw_cuda, ...
                        rot_angles_cuda, RotDet_cuda, param_cuda);
tmp = gather(Output_cuda);
gpuDevice(1); % reset gpu device
Output = zeros(size(q_all,1),9);
for j=1:size(q_all,1)
    Output(j,1:4) = tmp((j-1)*4+1:(j-1)*4+4);
    Output(j,5:9) = q_all(j,:);
end

Output_max=zeros(length(index_seeds(:,1)),9);
nr_per_OR1=OR_local.len;
Output=sortrows(Output,9,'ascend');
for j=1:length(index_seeds(:,1))
    if ~isempty(find(Output(:,9)==j))
%        ind0=find(Output(:,3)==max(Output(Output(:,9)==j,3)) & Output(:,9)==j);
        [~,ind0]=max(Output(1+nr_per_OR1*(j-1):nr_per_OR1*j,3));
        ind0=ind0+nr_per_OR1*(j-1);
        Output_max(j,:)=Output(ind0(1),:);
%         else
%             sprintf('Warning: OR %d returns empty result',j)
    end
end
Output_max = sortrows(Output_max,3,'descend');

