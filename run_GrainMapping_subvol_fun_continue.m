% continue running jobs for reconstruction subvolumes
% Jan 3, 2022
function run_GrainMapping_subvol_fun_continue(subvolnumber)

SampleFlag=3;
computation_options = [{'cpu'}, {'gpu_cuda_Gt'}, {'gpu_cuda_comp'}, {'gpu_cuda_index_compete'}];
compute_opt = computation_options{2};

if SampleFlag==2
    SampleName='AlCu8wt_middle_thinned_0930';
    % RecVolumePixel=[1   395;
    %                 1   399;
    %                 1   147]; % full volume
    RecVolumePixel_FOV=[63   311;
                    63   307;
                    45   115]; % effective volume masked by tomo
    OutputFolder=['./Examples/AlCu_8wt_middle_thinned_0930_rec/subvol_' num2str(subvolnumber)];
    % define total jobs = n(i)^3
    n=[3 3 1];
elseif SampleFlag==3
    SampleName='simu_Fe';
    RecVolumePixel_FOV=[1   160;
                1   160;
                1   240]; % effective volume masked by tomo
    OutputFolder=['./Examples/Fe_100um_11_11_simu_rec/subvol_' num2str(subvolnumber)];
    % define total jobs = n(i)^3
    n=[2 2 2];
elseif SampleFlag==4
    SampleName='virtual_Fe_100um_6grains';
    RecVolumePixel_FOV=[1   60;
                1   60;
                1   80]; % effective volume masked by tomo
    OutputFolder=['./Examples/virtual_Fe_100um_6grains_rec/subvol_' num2str(subvolnumber)];
    % define total jobs = n(i)^3
    n=[1 1 1];
else
    error('No such sample exists')
end

fname_prefix=['subvol_' num2str(subvolnumber)];
vol_divide_ind=[];
for i=1:3
    vol_divide_ind{i}=round(linspace(RecVolumePixel_FOV(i,1),RecVolumePixel_FOV(i,2),n(i)+1));
end

for i=1:n(1)
    for j=1:n(2)
        for k=1:n(3)
            ind=sub2ind([n(1) n(2) n(3)],i,j,k);
            RecVolumePixel_job{ind}(1,:)=vol_divide_ind{1}(i:i+1);
            RecVolumePixel_job{ind}(2,:)=vol_divide_ind{2}(j:j+1);
            RecVolumePixel_job{ind}(3,:)=vol_divide_ind{3}(k:k+1);
        end
    end
end
RecVolumePixel=RecVolumePixel_job{subvolnumber};
fprintf('The fullvolume has dimensions in X from %d to %d ...\n',RecVolumePixel(1,1),RecVolumePixel(1,2))
fprintf('The fullvolume has dimensions in Y from %d to %d ...\n',RecVolumePixel(2,1),RecVolumePixel(2,2))
fprintf('The fullvolume has dimensions in Z from %d to %d ...\n',RecVolumePixel(3,1),RecVolumePixel(3,2))
fprintf('Launching grain mapping for subvolume %d ...\n',subvolnumber)

if strcmp(compute_opt,'gpu_cuda_index_compete') == 1
    GrainMapping_index_and_compete_continue(OutputFolder,fname_prefix,RecVolumePixel,SampleName);
else
    GrainMapping_fun_continue(OutputFolder,fname_prefix,RecVolumePixel,SampleName,compute_opt);
end

