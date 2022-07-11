% run jobs for grain mapping as a whole volume
% Dec 6, 2021
% SampleName is consisten with setup_para.m
% e.g.: SampleName='AlCu8wt_middle_thinned_0930'; % 'simu_Fe'; 'AlCu8wt_middle_thinned_0930'
% computation_options: cpu, gpu_cuda_Gt, gpu_cuda_comp, gpu_cuda_index_compete
% update on April 8, 2022
function run_GrainMapping_fullvol_fun_continue(SampleFlag)
% SampleFlag=3;
computation_options = [{'cpu'}, {'gpu_cuda_Gt'}, {'gpu_cuda_comp'},{'gpu_cuda_index_compete'}];
compute_opt = computation_options{2};

if SampleFlag==2
    SampleName='AlCu8wt_middle_thinned_0930';
    % RecVolumePixel=[1   395;
    %                 1   399;
    %                 1   147]; % full volume
    RecVolumePixel_FOV=[63   311;
                    63   307;
                    45   115]; % effective volume masked by tomo
    OutputFolder=['./Examples/AlCu_8wt_middle_thinned_0930_rec/fullvol_' compute_opt];
    fname_prefix='fullvol';
elseif SampleFlag==3
    SampleName='simu_Fe';
    RecVolumePixel_FOV=[1   160;
                1   160;
                1   240]; % effective volume masked by tomo
   OutputFolder=['./Examples/Fe_100um_11_11_simu_rec/fullvol_' compute_opt];
   fname_prefix='fullvol';
elseif SampleFlag==4
    SampleName='virtual_Fe_100um_6grains';
    RecVolumePixel_FOV=[1   60;
                1   60;
                1   80]; % effective volume masked by tomo
   OutputFolder=['./Examples/virtual_Fe_100um_6grains_rec/test_' compute_opt];
   fname_prefix='fullvol';
else
    error('No such sample exists')
end
RecVolumePixel=RecVolumePixel_FOV;

sprintf('Continue grain reconstruction for %s ...',SampleName)
sprintf('The fullvolume has dimensions in X from %d to %d ...',RecVolumePixel(1,1),RecVolumePixel(1,2))
sprintf('The fullvolume has dimensions in Y from %d to %d ...',RecVolumePixel(2,1),RecVolumePixel(2,2))
sprintf('The fullvolume has dimensions in Z from %d to %d ...',RecVolumePixel(3,1),RecVolumePixel(3,2))

if strcmp(compute_opt,'gpu_cuda_index_compete') == 1
    GrainMapping_index_and_compete_continue(OutputFolder,fname_prefix,RecVolumePixel,SampleName);
else
    GrainMapping_fun_continue(OutputFolder,fname_prefix,RecVolumePixel,SampleName,compute_opt);
end


